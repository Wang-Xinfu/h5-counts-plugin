# dsh-tool-h5-counts

DeepSeek Harness 模型侧工具插件：读取 10x Genomics `.h5` 矩阵文件的**基因数 / 细胞数**（等价于 Harness 内置 `get_cell_gene_counts`）。

核心设计：**只读 `matrix/shape` 头部（8 字节），不加载、不解压矩阵数据** —— O(1) 读取，毫秒级返回规模信息，让 Agent「先看规模再决定要不要读全量」。

## 两个版本

本仓库提供同一工具的两个实现，行为一致（输入 `file_path`，输出基因数/细胞数），区别在读取路径：

| | 挂载版 `lib/index.js` | 完整版 `complete/get_cell_gene_counts.plugin.js` |
|---|---|---|
| 读取方式 | `python` + **h5py**（`pip install h5py`） | h5py 优先，**无 h5py 时退回内嵌纯标准库 HDF5 解析器** |
| 代码量 | ~70 行 | 721 行（内嵌 Python 读取器） |
| 兼容布局 | h5py 支持的一切 | superblock v0/v1/v2/v3、v1/v2 对象头、符号表 / Link 消息、chunked B-tree、deflate + shuffle 过滤器、HDF5 1.14 布局 |
| 依赖 | 本机 python + h5py | 仅 python（stdlib），h5py 可选加速 |
| 工程集成 | `execFile` 调 python | 脚本经 stdin 流入 `python -`、`H5_PATH` 环境变量传路径、fs 解析 / 沙箱 / 超时 / 中止 / 统一错误契约 |

## 安装（挂载版，挂到 web profile 宿主组合）

1. 克隆本仓库为 `C:\Users\MSI\.dsh\profiles\web\plugins\dsh-tool-h5-counts`：

   ```powershell
   git clone https://github.com/Wang-Xinfu/h5-counts-plugin.git C:\Users\MSI\.dsh\profiles\web\plugins\dsh-tool-h5-counts
   ```

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

## 完整版说明

`complete/get_cell_gene_counts.plugin.js` 是完整版宿主半区（`harness.registerTool` 契约，`inject: ['shell']`），内嵌 `README_SCRIPT`：

- **h5py 可用时优先**（覆盖所有布局）；
- **否则退回纯标准库读取器**：解析 superblock v0/v1（符号表）与 v2/v3（Link 消息，`OHDR` 魔数探测两种偏移排布）、v1/v2 对象头（含续接消息）、chunked 数据集的 B-tree（经典与 HDF5 1.14 布局逐节点校验）、deflate / shuffle 过滤器反向执行；对 `shape` 只读第一个 chunk 的前 64 字节即停；
- 文件路径经 `H5_PATH` 环境变量传递，脚本走 stdin 流入 `python -`，无 shell 转义问题；
- 每个字段做合法性校验（魔数、地址范围、rank 上限），失败输出 `{"error": "..."}` 并非零退出，宿主透传可操作提示（如 `run 'pip install h5py' and retry`）。

> 该文件按原样归档自开发工作区（`get_cell_gene_counts.plugin.js`），未做改动。若需接入不同版本的 Harness，按对应插件契约调整 `inject` 与注册方式即可。

## 依赖

- 挂载版运行期调用本机 `python` + `h5py`（`pip install h5py`）；完整版仅需 `python`（stdlib）。
- 若插件加载报 `tool "get_cell_gene_counts" is already registered in this scope`（与 Harness 内置工具重名），把 `lib/index.js` 中 `defineTool` 的 `name` 改成 `h5_matrix_counts` 即可。

## 性能基准

只读头部 vs 读取全量（本机实测，Python 3.14 + h5py 3.16，10x 官方 PBMC 3k 数据集
`pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix.h5`，38.8 MB，134,920 × 2,711，24,511,186 个非零元素）：

| 做法 | 读取量 | 耗时 |
|---|---|---|
| 只读 `matrix/shape`（本插件） | 8 字节 | **35 ms** |
| 读取完整矩阵（`data`+`indices`+`indptr`） | 2451 万个非零元素 | **520 ms** |

工作区构造的 v0/v2/chunked/shuffle 各布局测试矩阵，头部读取在 **0.77–2.12 ms**。

## 卸载

- 从 `cordis.patch.yml` 移除 insert 行，重启；删除 profile 依赖与 plugins 目录。
