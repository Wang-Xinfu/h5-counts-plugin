# dsh-tool-h5-counts

模型侧工具插件：读取 10x Genomics .h5 矩阵文件的基因数 / 细胞数（等价于 Harness 内置
`get_cell_gene_counts`，仅读 `matrix/shape` 头部，不解压矩阵数据）。

## 安装（挂到 web profile 宿主组合）

1. 把本目录复制为 `C:\Users\MSI\.dsh\profiles\web\plugins\dsh-tool-h5-counts`。
2. 编辑 `C:\Users\MSI\.dsh\profiles\web\package.json`，在 `dependencies` 中加入：

   ```json
   "dependencies": {
     "dsh-tool-h5-counts": "file:./plugins/dsh-tool-h5-counts"
   }
   ```

3. 在 profile 目录安装依赖（转发 pnpm）：

   ```powershell
   dsh plugin --profile web install
   ```

4. 编辑 `C:\Users\MSI\.dsh\profiles\web\cordis.patch.yml`，加入：

   ```yaml
   - insert:
       - id: tool-h5-counts
         name: dsh-tool-h5-counts
   ```

5. 重启 `dsh web`，新会话即可使用 `get_cell_gene_counts`。

## 依赖

- 运行期调用本机 `python` + `h5py`（`pip install h5py`）。
- 若插件加载报 `tool "get_cell_gene_counts" is already registered in this scope`，
  把 `lib/index.js` 中 `defineTool` 的 `name` 改成 `h5_matrix_counts` 即可。

## 卸载

- 从 `cordis.patch.yml` 移除 insert 行，重启；删除 profile 依赖与 plugins 目录。
