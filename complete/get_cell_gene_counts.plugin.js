/* =====================================================================
 * get_cell_gene_counts - Host half of the dynamic Cordis Plugin.
 *
 * Registers one model tool, get_cell_gene_counts:
 *   input : file_path - path to a 10x Genomics .h5 matrix file
 *   output: { gene_count, cell_count } (rows / columns of the matrix)
 *
 * The .h5 file is parsed by an embedded, dependency-free Python script
 * (h5_shape_reader.py) that reads only the matrix/shape dataset. The
 * script is streamed to `python -` over stdin; the file path travels in
 * the H5_PATH environment variable so no shell quoting is ever involved.
 *
 * h5py is used automatically when it is installed (handles every layout,
 * including chunked and compressed storage). The built-in stdlib reader
 * covers the layouts h5py / Cell Ranger actually produce: superblock
 * v0/v1 with v1 object headers (continuation-first, symbol-table message
 * 0x0011) and symbol tables; superblock v2/v3 with v2 object headers and
 * Link messages; and chunked datasets with the deflate and shuffle
 * filters (zlib), in both the classic and HDF5 1.14 layout encodings.
 * ===================================================================== */
const READER_SCRIPT = `#!/usr/bin/env python3
"""Read the gene x cell shape of a 10x Genomics .h5 matrix file.

Prefers h5py when it is importable (covers every layout). Falls back to a
small stdlib-only HDF5 reader that understands the layouts h5py / Cell Ranger
actually produce: superblock v0/v1 with v1 object headers and symbol tables,
and best-effort superblock v2/v3 with v2 object headers and Link messages.

Only the dataset 'matrix/shape' (a 2 x int32 array [n_genes, n_cells]) is
read; no matrix data is loaded.

The file path is taken from the H5_PATH environment variable (avoids any
shell-quoting issues with the caller):

    H5_PATH=filtered_feature_bc_matrix.h5 python h5_shape_reader.py

Prints one JSON object on stdout: {"gene_count": N, "cell_count": M}
On failure prints {"error": "..."} and exits nonzero.
"""

import json
import os
import struct
import sys
import zlib

try:
    import h5py  # noqa: F401  (probe only; imported again inside main)
    HAVE_H5PY = True
except Exception:
    HAVE_H5PY = False
# Setting H5_FORCE_STDLIB=1 skips h5py and uses the built-in reader (debug).
if os.environ.get('H5_FORCE_STDLIB'):
    HAVE_H5PY = False

_NULL = 0xFFFFFFFFFFFFFFFF


def _u16(b, o):
    return struct.unpack_from('<H', b, o)[0]


def _u32(b, o):
    return struct.unpack_from('<I', b, o)[0]


def _u64(b, o):
    return struct.unpack_from('<Q', b, o)[0]


def _u(b, o):
    """Little-endian unsigned of 8/4/2 bytes, whichever fits."""
    n = len(b) - o
    if n >= 8:
        return struct.unpack_from('<Q', b, o)[0]
    if n >= 4:
        return struct.unpack_from('<I', b, o)[0]
    return struct.unpack_from('<H', b, o)[0]


def _unshuffle(raw, esize):
    """Undo the HDF5 shuffle filter: within each block, bytes are grouped by
    position across elements (all first bytes, all second bytes, ...)."""
    if esize <= 1 or len(raw) < esize:
        return raw
    n = len(raw) // esize
    out = bytearray(len(raw))
    for i in range(esize):
        for j in range(n):
            out[j * esize + i] = raw[i * n + j]
    return bytes(out)


class H5ShapeReader:
    """Minimal HDF5 reader that locates 'matrix/shape' and reads its data."""

    def __init__(self, path):
        self.path = path
        self.f = open(path, 'rb')
        self.f_size = os.fstat(self.f.fileno()).st_size
        self.off_size = 8
        self.len_size = 8
        self.base = 0
        self.root_addr = 0
        self.oh_version = 1
        try:
            self._parse_superblock()
        except Exception:
            self.f.close()
            raise

    def close(self):
        self.f.close()

    def _read(self, addr, size):
        if addr in (0, _NULL):
            return b''
        self.f.seek(self.base + addr)
        return self.f.read(size)

    def _parse_superblock(self):
        self.f.seek(0)
        sig = self.f.read(8)
        if sig != b'\\x89HDF\\r\\n\\x1a\\n':
            raise ValueError('not an HDF5 file (bad signature)')
        ver = self.f.read(1)[0]
        if ver in (0, 1):
            self.f.seek(13)
            self.off_size = self.f.read(1)[0]
            self.f.seek(14)
            self.len_size = self.f.read(1)[0]
            if self.off_size not in (2, 4, 8) or self.len_size not in (2, 4, 8):
                raise ValueError(
                    'unsupported HDF5 address sizes %d/%d' % (self.off_size, self.len_size))
            self.f.seek(24)
            self.base = _u(self.f.read(self.off_size), 0)
            self.f.seek(24 + 4 * self.off_size + 8)  # root sym entry, oh addr at +8
            self.root_addr = _u(self.f.read(self.off_size), 0)
            self.oh_version = 1
        elif ver in (2, 3):
            self.off_size = 8
            self.len_size = 8
            self.base = 0
            # Modern v2/v3 superblocks (HDF5 1.10+) have no reserved byte:
            # base @12, root group OH @36. Some writers (and the v2 spec)
            # use a reserved byte: base @13, root @37. Validate both against
            # the OHDR signature and pick the layout that matches the file.
            self.f.seek(12)
            base_a = _u64(self.f.read(8), 0)
            self.f.seek(36)
            root_a = _u64(self.f.read(8), 0)
            self.f.seek(13)
            base_b = _u64(self.f.read(8), 0)
            self.f.seek(37)
            root_b = _u64(self.f.read(8), 0)

            def looks_like_ohdr(addr):
                if not 0 < addr < self.f_size:
                    return False
                self.f.seek(addr)
                return self.f.read(4) == b'OHDR'

            if looks_like_ohdr(root_a):
                self.base, self.root_addr = base_a, root_a
            elif looks_like_ohdr(root_b):
                self.base, self.root_addr = base_b, root_b
            else:
                self.base, self.root_addr = base_a, root_a
            self.oh_version = 2
        else:
            raise ValueError('unsupported HDF5 superblock version %d' % ver)

    # ---- public -----------------------------------------------------------

    def matrix_shape(self):
        """Return (gene_count, cell_count) from 'matrix/shape'."""
        members = self._group_members(self.root_addr)
        matrix = self._lookup(members, 'matrix')
        if matrix is None:
            raise ValueError(
                "could not find group 'matrix' in the file (not a 10x matrix .h5?)")
        members = self._group_members(matrix)
        shape = self._lookup(members, 'shape')
        if shape is None:
            raise ValueError("could not find dataset 'matrix/shape' in the file")
        dims = None
        layout = None
        filters = None
        for mtype, data in self._object_messages(shape):
            if mtype == 0x0001:
                dims = self._dataspace_dims(data)
            elif mtype == 0x0008:
                layout = data
            elif mtype == 0x000B:
                filters = data  # filter pipeline (needed for chunked data)
        if dims is None:
            raise ValueError("dataset 'matrix/shape' has no dataspace message")
        if layout is None:
            raise ValueError("dataset 'matrix/shape' has no data layout message")
        raw = self._read_layout(layout, dims, filters)
        if len(raw) < 8:
            raise ValueError('matrix/shape data is shorter than 8 bytes')
        n_genes, n_cells = struct.unpack_from('<ii', raw, 0)
        if n_genes < 0 or n_cells < 0:
            raise ValueError('matrix/shape contains negative dimensions')
        return n_genes, n_cells

    # ---- object headers ---------------------------------------------------

    def _object_messages(self, addr):
        if self.oh_version == 1:
            return self._v1_messages(addr)
        return self._v2_messages(addr)

    def _v1_messages(self, addr):
        hdr = self._read(addr, 16)
        if len(hdr) < 16 or hdr[0] != 1:
            raise ValueError('bad v1 object header at %d' % addr)
        hsize = _u32(hdr, 8)
        regions = [(addr + 16, hsize)]
        msgs = []
        while regions:
            start, size = regions.pop()
            pos = 0
            guard = 0
            while pos + 8 <= size and guard < 256:
                h = self._read(start + pos, 8)
                if len(h) < 8:
                    break
                mtype = _u16(h, 0)
                msize = _u16(h, 2)
                data = self._read(start + pos + 8, msize)
                if mtype == 0x0010:  # object header continuation (v1 layout)
                    if len(data) >= self.off_size + self.len_size:
                        coff = _u(data, 0)
                        clen = _u(data, self.off_size)
                        regions.append((coff, clen))
                else:
                    msgs.append((mtype, data))
                pos += 8 + ((msize + 7) & ~7)
                guard += 1
        return msgs

    def _v2_messages(self, addr):
        hdr = self._read(addr, 8)
        if len(hdr) < 8 or hdr[0:4] != b'OHDR' or hdr[4] != 2:
            raise ValueError('bad v2 object header at %d' % addr)
        flags = hdr[6]
        pos = 7  # after signature(4) + version(1) + reserved(1) + flags(1)
        # Optional fields selected by the header flags.
        if flags & 0x02:
            pos += 16   # creation/mod/access/change times
        if flags & 0x04:
            pos += 8    # attribute storage phase change values
        if flags & 0x08:
            pos += 8    # max compact / min dense attributes
        if flags & 0x10:
            pos += 4    # attribute dense flag
        pos += 2        # skip the 2-byte message count
        msgs = []
        stack = [addr + pos]
        while stack:
            start = stack.pop()
            pos = 0
            guard = 0
            in_block = 0
            while guard < 256:
                h = self._read(start + pos, 8)
                if len(h) < 8:
                    break
                mtype = _u16(h, 0)
                msize = _u16(h, 2)
                if mtype == 0 and msize == 0:
                    break  # Nil filler
                data = self._read(start + pos + 8, msize)
                if mtype == 0x0010 and len(data) >= 8:  # object header continuation
                    stack.append(_u64(data, 0))
                else:
                    msgs.append((mtype, data))
                pos += 8 + msize
                in_block += 1
                if in_block % 4 == 0:
                    pos += 4  # 4-byte gap after every four messages
                guard += 1
        return msgs

    # ---- groups (v1: symbol tables; v2: link messages) --------------------

    def _group_members(self, oh_addr):
        if self.oh_version == 1:
            return self._v1_group_members(oh_addr)
        return self._v2_group_members(oh_addr)

    def _v1_group_members(self, oh_addr):
        st = None
        for mtype, data in self._v1_messages(oh_addr):
            if mtype == 0x0011:  # symbol table message (v1 layout)
                st = data
                break
        if st is None:
            raise ValueError('group object header has no symbol table message')
        if len(st) < 2 * self.off_size:
            raise ValueError('bad symbol table message')
        btree_addr = _u(st, 0)
        heap_addr = _u(st, self.off_size)
        return self._btree_symbols(btree_addr, heap_addr)

    def _btree_symbols(self, btree_addr, heap_addr):
        out = []
        visited = set()
        snods = set()

        def walk(addr):
            if addr in visited or addr in (0, _NULL):
                return
            visited.add(addr)
            node = self._read(addr, 8 + 2 * self.off_size)
            if len(node) < 8 + 2 * self.off_size or node[0:4] != b'TREE':
                raise ValueError('bad group b-tree node at %d' % addr)
            level = node[5]
            nused = _u16(node, 6)
            if level == 0:
                esize = 40 + self.off_size
                for i in range(nused):
                    entry = self._read(
                        addr + 8 + 2 * self.off_size + i * esize, esize)
                    if len(entry) < esize:
                        continue
                    # The symbol-table-node address is stored in the leaf
                    # key's header field (entry+8) in real h5py files, in the
                    # key scratch pad (entry+24), and in the classic child
                    # pointer (entry+40). Every candidate is validated by the
                    # SNOD signature, so all layouts resolve correctly.
                    for cand in (_u(entry, 8), _u(entry, 24), _u(entry, 40)):
                        self._symbol_table_node(cand, heap_addr, out, snods)
            else:
                esize = 80 + self.off_size
                for i in range(nused):
                    entry = self._read(
                        addr + 8 + 2 * self.off_size + i * esize, esize)
                    if len(entry) < esize:
                        continue
                    walk(_u(entry, 80))

        walk(btree_addr)
        return out

    def _symbol_table_node(self, addr, heap_addr, out, visited):
        if addr in (0, _NULL) or addr in visited:
            return
        node = self._read(addr, 8)
        if len(node) < 8 or node[0:4] != b'SNOD':
            return  # not a symbol table node; not an error
        visited.add(addr)
        nsyms = _u16(node, 6)
        if nsyms > 4096:
            raise ValueError('suspicious symbol table node at %d' % addr)
        for i in range(nsyms):
            e = self._read(addr + 8 + i * 40, 40)
            if len(e) < 40:
                continue
            out.append((self._heap_name(heap_addr, _u64(e, 0)), _u(e, 8)))

    def _heap_name(self, heap_addr, name_off):
        hdr = self._read(heap_addr, 24 + self.off_size)
        if len(hdr) < 24 + self.off_size or hdr[0:4] != b'HEAP':
            raise ValueError('bad local heap at %d' % heap_addr)
        data_addr = _u(hdr, 24)
        if data_addr == 0:
            data_addr = heap_addr + 24 + self.off_size
        raw = self._read(data_addr + name_off, 4096)
        end = raw.find(b'\\x00')
        if end == -1:
            end = len(raw)
        return raw[:end].decode('utf-8', 'replace')

    def _v2_group_members(self, oh_addr):
        out = []
        try:
            messages = self._v2_messages(oh_addr)
        except Exception:
            raise ValueError(
                "the file uses an object-header variant the built-in reader "
                "cannot decode; run 'pip install h5py' and retry")
        for mtype, data in messages:
            if mtype == 0x0006:
                name, addr = self._link_target(data)
                if addr is not None:
                    out.append((name, addr))
            elif mtype == 0x0002:
                raise ValueError(
                    "group uses a fractal-heap link index the built-in reader "
                    "cannot decode; run 'pip install h5py' and retry")
        if not out:
            raise ValueError(
                "group object header has no link messages the built-in reader "
                "understands; run 'pip install h5py' and retry")
        return out

    @staticmethod
    def _link_target(data):
        if len(data) < 4:
            raise ValueError('bad link message')
        flags = data[1]
        name_len = _u16(data, 2)
        pos = 4
        if flags & 0x01:
            pos += 8    # creation order
        if flags & 0x02:
            pos += 1    # link type
        if flags & 0x04:
            pos += 1    # character set
        if pos + name_len > len(data):
            raise ValueError('bad link message name')
        name = data[pos:pos + name_len].decode('utf-8', 'replace')
        pos += name_len
        if pos + 8 > len(data):
            return name, None
        return name, _u64(data, pos)

    # ---- dataset messages -------------------------------------------------

    @staticmethod
    def _dataspace_dims(data):
        if len(data) < 4:
            raise ValueError('bad dataspace message')
        ver = data[0]
        rank = data[1]
        if rank > 32:
            raise ValueError('suspicious dataspace rank %d' % rank)
        if ver == 1:
            if len(data) < 4 + 8 * rank:
                raise ValueError('short dataspace message')
            return list(struct.unpack_from('<%dQ' % rank, data, 4))
        if ver == 2:
            if len(data) < 7 + 8 * rank:
                raise ValueError('short dataspace message')
            return list(struct.unpack_from('<%dQ' % rank, data, 7))
        raise ValueError('unsupported dataspace version %d' % ver)

    def _read_layout(self, data, dims, filters):
        if len(data) < 2:
            raise ValueError('bad data layout message')
        ver = data[0]
        cls = data[1]
        if cls == 0:  # compact: data lives inside the message
            if ver <= 2:
                size = _u16(data, 8)
                return data[10:10 + size]
            size = _u16(data, 2)
            return data[4:4 + size]
        if cls == 1:  # contiguous: address + size
            if ver <= 2:
                addr = _u(data, 8)
                size = _u(data, 8 + self.off_size)
            else:
                # HDF5 1.14 layout v3 has no reserved bytes: address at +2
                addr = _u64(data, 2)
                size = _u64(data, 10)
            return self._read(addr, min(size, 64))
        if cls == 2:  # chunked: read through the chunk b-tree, then filters
            return self._read_chunked(data, ver, dims, filters)
        raise ValueError('unsupported data layout class %d' % cls)

    # ---- chunked datasets -------------------------------------------------

    def _read_chunked(self, data, ver, dims, filters):
        rank = len(dims)
        if rank > 32:
            raise ValueError('suspicious chunked rank %d' % rank)
        if ver <= 2:
            if len(data) < 8 + 8 * rank + 8:
                raise ValueError('bad chunked layout message')
            btree_addr = _u(data, 8 + 8 * rank)
        else:
            # HDF5 1.14 layout v3 chunked: b-tree address at offset 3
            if len(data) < 11:
                raise ValueError('bad chunked layout message')
            btree_addr = _u64(data, 3)
        chunks = self._chunk_btree(btree_addr, rank, dims)
        if not chunks:
            raise ValueError('chunked b-tree holds no chunks')
        chunks.sort(key=lambda c: c[0])
        out = b''
        for _offs, csize, cmask, caddr in chunks:
            raw = self._read(caddr, csize)
            if len(raw) < csize:
                raise ValueError('chunk data truncated at %d' % caddr)
            out += self._apply_filters(raw, cmask, filters)
            if len(out) >= 64:
                break
        return out[:64]

    def _chunk_btree(self, btree_addr, rank, dims):
        """Collect (chunk_offset_tuple, size, filter_mask, address) for every
        leaf chunk of a chunked-dataset b-tree (node type 1).

        Handles both the classic entry layout (offset first, 4-byte size and
        filter mask, 20-28 bytes per entry) and the HDF5 1.14 layout (8-byte
        size first, element offset, 32 bytes per entry). The format is
        detected per node by validating the size/address fields."""
        idx_size = 8 if any(d > 0xFFFFFFFF for d in dims) else 4
        out = []
        visited = set()

        def entry_vals(entry):
            """Return ((offset_tuple, size, mask, addr), entry_size) or None."""
            # classic: offset(idx), size(4), mask(4), child(8)
            c_off = list(struct.unpack_from(
                '<%d%s' % (rank, 'I' if idx_size == 4 else 'Q'), entry, 0))
            c_size = _u32(entry, idx_size)
            c_mask = _u32(entry, idx_size + 4)
            c_addr = _u(entry, idx_size + 8)
            if 0 < c_size <= 64 * 1024 * 1024 and 0 < c_addr < self.f_size:
                return ((tuple(c_off), c_size, c_mask, c_addr), idx_size + 16)
            # HDF5 1.14: size(8), element offset(8), mask(4), pad(4), child(8)
            if len(entry) >= 32:
                n_size = _u64(entry, 0)
                n_addr = _u64(entry, 24)
                n_mask = _u32(entry, 16)
                if 0 < n_size <= 64 * 1024 * 1024 and 0 < n_addr < self.f_size:
                    return (((_u64(entry, 8),), n_size, n_mask, n_addr), 32)
            return None

        def walk(addr):
            if addr in visited or addr in (0, _NULL):
                return
            visited.add(addr)
            node = self._read(addr, 8 + 2 * self.off_size)
            if len(node) < 8 + 2 * self.off_size or node[0:4] != b'TREE':
                raise ValueError('bad chunked b-tree node at %d' % addr)
            level = node[5]
            nused = _u16(node, 6)
            if level == 0:
                base = addr + 8 + 2 * self.off_size
                probe = self._read(base, 40)
                got = entry_vals(probe) if len(probe) >= 16 else None
                if got is None:
                    raise ValueError('cannot decode chunk b-tree leaf at %d' % addr)
                _kv, esize = got
                for i in range(nused):
                    entry = self._read(base + i * esize, esize)
                    if len(entry) < esize:
                        continue
                    kv = entry_vals(entry)
                    if kv is not None:
                        out.append(kv[0])
            else:
                # classic internal entries: key = two chunk offsets + size +
                # mask, then the child address; never needed for matrix/shape
                esize = 2 * idx_size + 8 + self.off_size
                for i in range(nused):
                    entry = self._read(
                        addr + 8 + 2 * self.off_size + i * esize, esize)
                    if len(entry) < esize:
                        continue
                    walk(_u(entry, 2 * idx_size + 8))

        walk(btree_addr)
        return out

    def _apply_filters(self, raw, cmask, filters):
        """Reverse the filter pipeline stored on the chunk (mask bits = filters
        that were NOT applied to this chunk)."""
        if filters is None or len(filters) < 8:
            if cmask:
                raise ValueError(
                    'chunk is filtered but the dataset has no filter pipeline')
            return raw
        nfilters = filters[1]
        pos = 8  # version(1) + nfilters(1) + reserved(2) + 4-byte pad
        applied = []
        for i in range(nfilters):
            if pos + 16 > len(filters):
                break
            fid = _u16(filters, pos)
            name_len = _u16(filters, pos + 2)
            flags = filters[pos + 4]
            nvalues = _u16(filters, pos + 6)
            values = filters[pos + 8 + name_len:pos + 8 + name_len + 4 * nvalues]
            pos += 8 + name_len + 4 * nvalues
            pos += (8 - pos % 8) % 8  # each filter is padded to 8 bytes
            applied.append((fid, flags, values))
        for i in range(len(applied) - 1, -1, -1):
            fid, _flags, values = applied[i]
            if cmask & (1 << i):
                continue  # this filter was skipped for this chunk
            if fid == 1:  # deflate (zlib stream)
                raw = zlib.decompress(raw)
            elif fid == 2:  # shuffle (reorder bytes within elements)
                esize = _u32(values, 0) if len(values) >= 4 else 4
                raw = _unshuffle(raw, esize)
            else:
                raise ValueError(
                    'chunk uses filter id %d which the built-in reader cannot '
                    "decode; run 'pip install h5py' and retry" % fid)
        return raw

    @staticmethod
    def _lookup(members, name):
        for n, addr in members:
            if n == name:
                return addr
        return None


def main():
    path = os.environ.get('H5_PATH', '')
    if not path:
        print(json.dumps({'error': 'H5_PATH environment variable must name the 10x .h5 file'}))
        return 2
    try:
        if HAVE_H5PY:
            with h5py.File(path, 'r') as f:
                arr = f['matrix/shape'][...]
                n_genes = int(arr[0])
                n_cells = int(arr[1])
        else:
            reader = H5ShapeReader(path)
            try:
                n_genes, n_cells = reader.matrix_shape()
            finally:
                reader.close()
        print(json.dumps({'gene_count': n_genes, 'cell_count': n_cells}))
        return 0
    except Exception as exc:
        print(json.dumps({'error': '%s: %s' % (type(exc).__name__, exc)}))
        return 1


if __name__ == '__main__':
    sys.exit(main())
`;

return {
  name: 'get-cell-gene-counts',
  inject: ['shell'],
  apply(ctx) {
    harness.registerTool(ctx, harness.defineTool({
      name: 'get_cell_gene_counts',
      description: 'Read a 10x Genomics .h5 single-cell matrix file and return its cell count (number of columns) and gene count (number of rows). The file must be an HDF5 (.h5) file in 10x matrix format (e.g. filtered_feature_bc_matrix.h5 from Cell Ranger or a Scanpy write_10x_h5 output). Chunked and compressed storage is supported (h5py when installed; otherwise the built-in reader handles deflate and shuffle filters). Only the stored matrix/shape header is read; the matrix data itself is never loaded.',
      parameters: {
        file_path: {
          type: 'string',
          required: true,
          description: 'Path to the 10x Genomics .h5 matrix file to inspect (absolute, or relative to the session workspace).',
        },
      },
      output: {
        schema: {
          type: 'object',
          additionalProperties: false,
          properties: {
            cell_count: { type: 'integer', required: true, description: 'Number of cells = number of columns in the expression matrix.' },
            gene_count: { type: 'integer', required: true, description: 'Number of genes = number of rows in the expression matrix.' },
          },
        },
        render: (_args, value) => [{
          type: 'text',
          text: `get_cell_gene_counts: ${value.gene_count} genes (rows) x ${value.cell_count} cells (columns)`,
        }],
      },
      async execute(args, exec) {
        if (typeof args.file_path !== 'string' || args.file_path.trim() === '') {
          throw new Error('invalid file_path: expected a non-empty string path to a 10x .h5 file')
        }
        const headerCwd = exec.agent ? exec.agent.session.header.cwd : undefined
        const fs = ctx.get('fs')
        let processPath = args.file_path
        if (fs !== undefined) {
          const target = await fs.resolve(args.file_path, headerCwd === undefined ? undefined : { cwd: headerCwd })
          const info = await fs.stat(target)
          if (info === undefined) throw new Error(`file not found: ${args.file_path}`)
          processPath = fs.processPath(target)
        }
        const request = {
          command: 'python -',
          stdin: READER_SCRIPT,
          timeoutMs: 60000,
          stdoutMaxBytes: 262144,
        }
        if (headerCwd !== undefined) request.workdir = headerCwd
        const shellEnv = ctx.get('shellEnv')
        if (shellEnv !== undefined) request.dshEnv = shellEnv.collect(exec)
        const sandboxPolicy = ctx.get('sandboxPolicy')
        if (ctx.shell.sandboxMode !== undefined && sandboxPolicy !== undefined) {
          request.sandboxPolicy = sandboxPolicy.resolve(exec.agent === undefined ? {} : { session: exec.agent.session })
        }
        request.env = { H5_PATH: processPath }
        const result = await ctx.shell.run(ctx.shell.resolve({ ...request, signal: exec.signal }))
        if (result.aborted) {
          const error = new Error('tool call aborted')
          error.name = 'AbortError'
          throw error
        }
        if (result.sandbox && result.sandbox.denied) {
          throw new Error(`[sandbox: file access denied under ${result.sandbox.mode} mode] the .h5 file is outside the permitted workspace: ${args.file_path}`)
        }
        let parsed
        try {
          parsed = JSON.parse(result.stdout.text)
        } catch (error) {
          const detail = result.stderr.text.trim()
          throw new Error(`get_cell_gene_counts: python reader failed (exit ${result.exitCode})${detail ? ': ' + detail : ''}`)
        }
        if (parsed && parsed.error) throw new Error(`get_cell_gene_counts: ${parsed.error}`)
        if (!parsed || !Number.isInteger(parsed.gene_count) || !Number.isInteger(parsed.cell_count)) {
          throw new Error(`get_cell_gene_counts: unexpected reader output ${JSON.stringify(parsed)}`)
        }
        return { cell_count: parsed.cell_count, gene_count: parsed.gene_count }
      },
      presentCall: (args) => ({
        card: 'generic',
        title: 'get_cell_gene_counts',
        kind: 'execute',
        rawInput: String(args.file_path),
        content: [{ type: 'text', text: `Reading 10x matrix: ${args.file_path}` }],
      }),
    }))
  },
}