# IQ-TREE Agent Skills

面向 LLM Agent 的 IQ-TREE 技能集，内容整理自官方文档 <https://iqtree.github.io/doc/>。
每个技能文件夹可直接装载进支持 Agent Skills 规范的 Agent 工具（如 Claude Code、Qoder 等）：
Agent 读取 `SKILL.md` 获得该主题的工作指引，再按需加载 `references/` 中的官方文档镜像。

## 技能列表（12 个）

| 技能文件夹 | 主题 | 参考文档数 |
| --- | --- | --- |
| `iqtree/` | 入口技能：安装、首次推断、输入格式与数据类型、多核、输出文件解读、checkpoint 续算、大数据策略、Web 服务器、FAQ 排错、源码编译 | 9 |
| `iqtree-alisim/` | 序列模拟（AliSim）：随机树/固定树模拟、Indel、异速演化（GHOST）、多数据类型混合 | 2 |
| `iqtree-branch-support/` | 分支支持度：UFBoot、非参数 bootstrap、SH-aLRT/aBayes/LBP、对称检验与似然映射 | 2 |
| `iqtree-command-reference/` | 命令行逐项参考：全部选项、1.x/2.x/3.x 拼写对照、输出文件清单，附 3 个辅助脚本 | 3 |
| `iqtree-complex-models/` | 复杂模型：分区模型（-p/-q/-Q）、合并方案、混合模型 MIX/FMIX、PMSF、GHOST、MAST | 2 |
| `iqtree-concordance-factor/` | 演化不一致分析：gCF/sCF、concordance vectors、序列放置（phylogenetic placement） | 3 |
| `iqtree-dating/` | 分子钟定年：IQ2MC/MCMCtree 工作流、LSD2 全部场景、`--date-options` 原生参数 | 2 |
| `iqtree-model-selection/` | 模型选择：ModelFinder、DNA/蛋白/密码子/二态/形态模型目录、QMaker/nQMaker 自定义矩阵 | 2 |
| `iqtree-poms/` | PoMo 多态性感知模型：counts 文件、+P/+N/+WB/+WH 语法、物种树推断 | 1 |
| `iqtree-recipes/` | 端到端实战配方：序列放置、Goldman-Cox 检验、concordance vectors、多树连接模拟 | 6 |
| `iqtree-rooting/` | 无外群定根：Rootstrap、非可逆模型（NONREV/NQ.*）、`--root-test` | 2 |
| `iqtree-tree-search/` | 树搜索控制：起始/固定树、约束搜索 -g、NNI/SPR 参数、病原体数据快速搜索、拓扑检验 | 2 |

`iqtree/` 是推荐入口：其 "Routing table" 一节按用户目标把任务路由到其余技能。

## 目录结构

```
.
├── iqtree/SKILL.md          # 每个技能一份入口指南（YAML frontmatter: name + description）
│   └── references/*.md      # 官方文档逐字镜像（保留原 frontmatter，含 source_url 溯源）
├── iqtree-command-reference/
│   └── scripts/             # build_command.py / check_env.py / parse_results.py（Python 3，无第三方依赖）
└── assets/
    ├── images/              # 官方文档图片（14 个）
    └── files/               # 文档引用的示例数据（bovidae、plant、astral 等，12 个）
```

## 使用方法

- **装载技能**：把所需技能文件夹复制到你的 Agent 技能目录（例如 Claude Code 的
  `.claude/skills/<skill-name>/`）。技能之间通过相对路径 `../<skill>/references/...`
  互相引用，建议整仓一起使用。
- **只取所需**：`references/` 是整页镜像，各 SKILL.md 的 "grep recipes" 一节给出按标题
  精确抽取片段的命令，避免整页加载。
- **辅助脚本**（`iqtree-command-reference/scripts/`，依赖仅 Python 3 标准库）：

```bash
python3 scripts/check_env.py                 # 找到 iqtree/iqtree2/iqtree3 二进制并推荐 -T/-mem
python3 scripts/build_command.py --alignment aln.phy -T 4 --bootstrap 1000 --alrt 1000
python3 scripts/parse_results.py --prefix run1
# 三个脚本均支持 --self-test；build_command.py 会校验对齐文件存在，
# 文件尚未生成时可加 --allow-missing 先产出命令行
```

## 数据来源与许可

- `references/` 下的文档镜像与 `assets/` 下的图片、示例数据来自官方文档站
  [iqtree.github.io/doc](https://iqtree.github.io/doc/)，版权归 IQ-TREE 作者所有，
  此处仅作本地镜像以便离线检索；每个文件 frontmatter 中的 `source_url` 记录了原始出处。
- `SKILL.md` 指南、脚本与仓库组织采用 [MIT 许可](LICENSE)。
- 引用 IQ-TREE 时请使用官方首页列出的文献（见 `iqtree/references/home.md`）。
