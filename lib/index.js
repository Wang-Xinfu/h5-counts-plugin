// dsh-tool-h5-counts —— 模型侧 10x .h5 计数工具插件
// 与 Harness 内置 get_cell_gene_counts 行为一致：
// 只读取 /matrix/shape 头部，不加载、不解压矩阵数据（毫秒级）。
// 读取复用本机 Python + h5py；若换机器需保证 python 与 h5py 可用。
import { execFile } from 'node:child_process';
import { promisify } from 'node:util';
import { defineTool } from '@deepseek-ai/dsh-tools';

const execFileAsync = promisify(execFile);

// Cordis 插件契约
const name = 'tool-h5-counts';
const inject = ['tools'];

const PY_SCRIPT = [
  'import json, sys',
  'import h5py',
  "with h5py.File(sys.argv[1], 'r') as f:",
  "    shape = f['matrix/shape'][()]",
  "print(json.dumps({'genes': int(shape[0]), 'cells': int(shape[1])}))",
].join('\n');

async function runPython(path) {
  try {
    const { stdout } = await execFileAsync('python', ['-c', PY_SCRIPT, path], {
      maxBuffer: 4 * 1024 * 1024,
      windowsHide: true,
    });
    return JSON.parse(stdout.trim());
  } catch (error) {
    throw new Error(
      `cannot read "${path}" as a 10x .h5 matrix (python/h5py must be on PATH): ${error.message}`,
    );
  }
}

function apply(ctx) {
  ctx.tools.register(defineTool({
    name: 'get_cell_gene_counts',
    description:
      'Read a 10x Genomics .h5 single-cell matrix file and return its cell count (number of columns) and gene count (number of rows). The file must be an HDF5 (.h5) file in 10x matrix format (e.g. filtered_feature_bc_matrix.h5 from Cell Ranger or a Scanpy write_10x_h5 output). Only the stored matrix/shape header is read; the matrix data itself is never loaded.',
    parameters: {
      file_path: {
        type: 'string',
        required: true,
        description:
          'Path to the 10x Genomics .h5 matrix file to inspect (absolute, or relative to the session workspace).',
      },
    },
    output: {
      schema: {
        type: 'object',
        properties: {
          genes: { type: 'integer', required: true, description: 'Gene count (matrix rows).' },
          cells: { type: 'integer', required: true, description: 'Cell count (matrix columns).' },
        },
        additionalProperties: false,
      },
      render(_args, value) {
        return [{
          type: 'text',
          text: `get_cell_gene_counts: ${value.genes} genes (rows) x ${value.cells} cells (columns)`,
        }];
      },
    },
    async execute(args) {
      return runPython(args.file_path);
    },
  }));
}

export { name, inject, apply };
