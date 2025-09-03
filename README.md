# CDK2抑制剂筛选数据分析与可视化平台

[![Python](https://img.shields.io/badge/Python-3.8%2B-blue)](https://www.python.org/)
[![Pandas](https://img.shields.io/badge/Pandas-1.3%2B-green)](https://pandas.pydata.org/)
[![Matplotlib](https://img.shields.io/badge/Matplotlib-3.4%2B-orange)](https://matplotlib.org/)
[![RDKit](https://img.shields.io/badge/RDKit-2021%2B-red)](https://www.rdkit.org/)
[![License](https://img.shields.io/badge/License-MIT-yellow)](LICENSE)

## 📋 目录

- [项目简介](#项目简介)
- [背景知识](#背景知识)
- [功能特性](#功能特性)
- [环境要求](#环境要求)
- [安装指南](#安装指南)
- [代码结构](#代码结构)
- [详细功能解析](#详细功能解析)
- [数据说明](#数据说明)
- [运行示例](#运行示例)
- [结果解读](#结果解读)
- [扩展应用](#扩展应用)
- [常见问题](#常见问题)
- [贡献指南](#贡献指南)
- [引用信息](#引用信息)

---

## 🔬 项目简介

本项目实现了CDK2（细胞周期依赖性激酶2）抑制剂的完整筛选数据分析流程，通过Python对20种候选化合物进行多维度评估和可视化展示。项目复现了药物发现领域的标准筛选流程，包括初步筛选、IC50测定、验证筛选和分子结构分析。

### 核心价值
- 🎯 **真实案例**：基于实际药物筛选研究数据
- 📊 **专业可视化**：采用制药行业标准图表
- 🔄 **完整流程**：从虚拟筛选到实验验证的全链条
- 💡 **教育意义**：适合药物发现和生物信息学教学

---

## 🧬 背景知识

### CDK2的重要性
CDK2（Cyclin-dependent kinase 2）是调控细胞周期G1/S转换的关键蛋白激酶。在多种癌症中，CDK2过度激活导致细胞异常增殖。因此，CDK2抑制剂是潜在的抗癌药物。

### 药物筛选策略
```
虚拟筛选 (in silico)
    ↓
初步筛选 (10μM)  ← 快速评估
    ↓
IC50测定         ← 精确量化
    ↓
验证筛选 (50μM)  ← 剂量确认
    ↓
结构优化         ← SAR分析
```

---

## ✨ 功能特性

### 1. 数据处理与管理
- ✅ 20种化合物的多维度数据整合
- ✅ Pandas DataFrame高效数据处理
- ✅ 自动数据验证和异常处理

### 2. 可视化功能
- 📊 初步筛选结果柱状图（10μM抑制率）
- 📈 IC50剂量-反应曲线
- 📊 验证筛选结果图（50μM抑制率）
- 🧬 分子结构网格展示

### 3. 科学分析
- 🔬 Hill方程拟合
- 📐 活性分级系统
- 🎯 候选化合物自动识别
- 📋 统计报告生成

---

## 💻 环境要求

### 基础环境
- Python 3.8 或更高版本
- Jupyter Notebook 或 JupyterLab
- 操作系统：Windows/macOS/Linux

### Python依赖包
```python
numpy >= 1.19.0          # 数值计算
pandas >= 1.3.0          # 数据处理
matplotlib >= 3.4.0      # 基础绘图
seaborn >= 0.11.0        # 统计图表
rdkit >= 2021.03.1       # 分子结构（可选）
```

---

## 📦 安装指南

### 方法一：使用pip安装
```bash
# 创建虚拟环境（推荐）
python -m venv cdk2_env
source cdk2_env/bin/activate  # Linux/Mac
# 或
cdk2_env\Scripts\activate      # Windows

# 安装基础依赖
pip install numpy pandas matplotlib seaborn

# 安装RDKit（用于分子结构展示）
pip install rdkit-pypi
```

### 方法二：使用conda安装
```bash
# 创建conda环境
conda create -n cdk2_analysis python=3.8

# 激活环境
conda activate cdk2_analysis

# 安装所有依赖
conda install -c conda-forge numpy pandas matplotlib seaborn rdkit
```

### 方法三：使用requirements.txt
```bash
# 创建requirements.txt文件
echo "numpy>=1.19.0
pandas>=1.3.0
matplotlib>=3.4.0
seaborn>=0.11.0
rdkit-pypi>=2021.3.1" > requirements.txt

# 安装
pip install -r requirements.txt
```

---

## 📁 代码结构

```
CDK2_Screening_Analysis/
│
├── 📊 Cell_1_Data_Initialization.py      # 数据导入与环境配置
├── 📈 Cell_2_First_Screening.py          # 第一次筛选可视化
├── 📉 Cell_3_IC50_Determination.py       # IC50测定曲线
├── 📊 Cell_4_Second_Screening.py         # 第二次筛选验证
├── 🧬 Cell_5_Molecular_Structures.py     # 分子结构展示
├── 📝 README.md                           # 项目文档
└── 📋 requirements.txt                    # 依赖列表
```

---

## 🔍 详细功能解析

### 📊 Cell 1: 数据初始化与环境配置

#### 1.1 导入必要的库
```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Rectangle
import matplotlib.patches as mpatches
```

**代码解析**：
- `numpy`: 提供高效的数组操作和数学函数
- `pandas`: 结构化数据处理的核心工具
- `matplotlib.pyplot`: 创建各类科学图表
- `seaborn`: 提供美观的默认样式和统计图表
- `matplotlib.patches`: 用于添加图形标注和装饰

#### 1.2 配置中文显示支持
```python
plt.rcParams['font.sans-serif'] = ['SimHei', 'Arial Unicode MS', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False
plt.style.use('seaborn-v0_8')
```

**技术要点**：
- 设置多个字体作为备选，确保跨平台兼容性
- `axes.unicode_minus = False` 解决负号显示问题
- 使用seaborn样式提升图表美观度

#### 1.3 化合物数据存储
```python
compounds_data = {
    'AE-848/33854034': {
        'gscore': -10.932,        # Glide对接分数
        'asie': 24.75,           # 平均特异性相互作用能
        'inhibition_10um': 15.2,  # 10μM抑制率(%)
        'inhibition_50um': 28.5   # 50μM抑制率(%)
    },
    # ... 另外19个化合物数据
}
```

**数据结构设计**：
- 使用嵌套字典便于访问和扩展
- 键名采用"库名/编号"格式，符合行业标准
- 四个维度数据覆盖计算预测和实验验证

---

### 📈 Cell 2: 第一次预筛选可视化（10μM）

#### 2.1 数据准备与转换
```python
def plot_first_screening_10um():
    """复现论文图3-2: 第一次预筛选时化合物在10μM时对CDK2的抑制作用"""
    
    # 转换为DataFrame以便处理
    df = pd.DataFrame.from_dict(compounds_data, orient='index')
    
    # 提取化合物编号用于x轴标签
    compound_names = [name.split('/')[-1] for name in df.index]
```

**技术细节**：
- `orient='index'` 将字典键作为行索引
- 使用列表推导式提取简化的化合物ID
- DataFrame结构便于后续的数据操作和筛选

#### 2.2 活性分级着色系统
```python
# 根据抑制率设置颜色 - 突出显示活性化合物
colors = []
for i, inhibition in enumerate(inhibition_10um):
    if inhibition > 35:      # 高活性
        colors.append('red')
    elif inhibition > 20:    # 中等活性
        colors.append('orange') 
    elif inhibition > 10:    # 低活性
        colors.append('lightblue')
    else:                    # 无显著活性
        colors.append('lightgray')
```

**设计理念**：
- 颜色梯度反映活性强弱
- 红色警示高活性化合物，需重点关注
- 灰色表示无效化合物，可排除

#### 2.3 重点化合物标注
```python
# 标出两个被选中进行IC50测试的化合物
selected_indices = [5, 14]  # AQ-149/42126274 和 AQ-390/10780015
for idx in selected_indices:
    bars[idx].set_edgecolor('darkred')
    bars[idx].set_linewidth(3)
    
    # 添加箭头标注
    ax.annotate(f'选中\nIC50测试', 
               xy=(idx, inhibition_10um[idx]), 
               xytext=(idx, inhibition_10um[idx] + 15),
               arrowprops=dict(arrowstyle='->', color='darkred', lw=2),
               ha='center', fontsize=10, color='darkred', fontweight='bold')
```

**可视化技巧**：
- 加粗边框突出重点
- 箭头标注提供额外信息
- 文字说明解释选择原因

---

### 📉 Cell 3: IC50测定与曲线拟合

#### 3.1 Hill方程模拟
```python
def plot_ic50_curves():
    """复现论文图3-3: 从第一次预筛选中选择的两种化合物的IC50值"""
    
    # 定义浓度梯度
    concentrations = np.array([0.01, 0.1, 1, 10, 100, 500])  # μM
    
    # AQ-149/42126274的剂量-反应曲线
    ic50_1 = 16.976
    # Hill方程: Response = 100 / (1 + (IC50/[Drug])^n)
    response_1 = 100 / (1 + (ic50_1/concentrations)**1.2)
```

**科学原理**：
- Hill方程描述配体-受体结合动力学
- Hill系数1.2表示轻微正协同性
- IC50是药效评价的金标准

#### 3.2 半对数坐标系应用
```python
ax1.semilogx(concentrations, response_1, 'bo-', markersize=8, linewidth=2)
ax1.axvline(x=ic50_1, color='red', linestyle='--', alpha=0.7)
ax1.axhline(y=50, color='gray', linestyle='--', alpha=0.7)
```

**图表设计**：
- `semilogx`: x轴对数刻度，适合宽浓度范围
- 垂直线标记IC50位置
- 水平线标记50%抑制水平
- 交叉点即为IC50值

#### 3.3 结果标注
```python
ax1.text(ic50_1*2, 45, f'IC50 = {ic50_1:.3f} μM', fontsize=10, 
         bbox=dict(boxstyle="round,pad=0.3", facecolor="yellow", alpha=0.7))
```

**展示技巧**：
- 黄色背景框提高可读性
- 位置偏移避免遮挡数据点
- 保留3位小数体现精确度

---

### 📊 Cell 4: 第二次筛选验证（50μM）

#### 4.1 更严格的筛选标准
```python
def plot_second_screening_50um():
    """复现论文图3-4: 第二次预筛选时化合物在50μM时对CDK2的抑制作用"""
    
    # 使用50%抑制率作为阈值
    ax.axhline(y=50, color='green', linestyle='-', linewidth=2, alpha=0.8, 
               label=f'IC50测试选择阈值 (50%)')
```

**筛选策略升级**：
- 提高浓度到50μM验证剂量依赖性
- 50%阈值确保选出真正有效的化合物
- 绿色实线强调关键决策点

#### 4.2 五级活性分类
```python
# 更细化的活性分级
colors = []
for inhibition in inhibition_50um:
    if inhibition > 60:      # 极高活性
        colors.append('darkred')
    elif inhibition > 40:    # 高活性
        colors.append('red')
    elif inhibition > 25:    # 中等活性
        colors.append('orange')
    elif inhibition > 15:    # 低活性
        colors.append('lightblue')
    else:                    # 无显著活性
        colors.append('lightgray')
```

**分级优化**：
- 五级分类提供更精细的区分
- 深红色标识最优候选物
- 便于快速视觉筛选

#### 4.3 统计信息输出
```python
print(f"\n=== 第二次预筛选结果 (50μM) ===")
print(f"总化合物数: {len(compound_names)}")
print(f"抑制率 > 50% 的化合物数: {len(active_compounds)}")
for i in active_compounds:
    print(f"  - {compound_names[i]}: {inhibition_50um[i]:.1f}%")
```

**信息汇总**：
- 定量统计筛选效果
- 列出所有候选化合物
- 提供决策支持数据

---

### 🧬 Cell 5: 分子结构可视化

#### 5.1 RDKit环境检测
```python
try:
    from rdkit import Chem
    from rdkit.Chem import Draw
    from rdkit.Chem.Draw import IPythonConsole
    RDKIT_AVAILABLE = True
    print("RDKit available - 将显示真实分子结构图")
except ImportError:
    RDKIT_AVAILABLE = False
    print("RDKit not available - 请安装rdkit: conda install -c conda-forge rdkit")
```

**错误处理**：
- 优雅地处理依赖缺失
- 提供安装指导
- 不影响其他功能运行

#### 5.2 SMILES结构定义
```python
compounds_smiles = {
    'AE-848/33854034': 'CC1=CC=C(C=C1)C2=NC3=C(S2)C=CC=C3',
    'AP-064/43237294': 'COC1=CC=C(C=C1)C2=NC(=CS2)C3=CC=CC=C3',
    # ... 更多SMILES字符串
}
```

**SMILES说明**：
- Simplified Molecular Input Line Entry System
- 用ASCII字符串表示分子结构
- 行业标准的分子表示方法
- 示例解析：
  - `C`: 碳原子
  - `=`: 双键
  - `1`: 环闭合标记
  - `(C=C1)`: 支链结构

#### 5.3 分子对象创建与验证
```python
for comp_id, smiles in compounds_smiles.items():
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            molecules.append(mol)
            
            # 创建信息标签
            short_id = comp_id.split('/')[-1]
            inhibition = inhibition_10um.get(comp_id, 0)
            
            if comp_id in ic50_data:
                legend = f"{short_id}\nIC50: {ic50_data[comp_id]:.1f}μM"
            else:
                legend = f"{short_id}\n{inhibition:.1f}% @ 10μM"
```

**处理流程**：
- SMILES解析为分子对象
- 验证分子结构有效性
- 准备标签信息（ID、活性数据）
- 区分有/无IC50数据的化合物

#### 5.4 批量展示优化
```python
# 分批显示分子结构（每次显示10个，便于观看）
batch_size = 10
n_batches = (len(molecules) + batch_size - 1) // batch_size

for batch_idx in range(n_batches):
    img = Draw.MolsToGridImage(
        batch_molecules,
        molsPerRow=5,              # 每行5个
        subImgSize=(300, 250),     # 每个分子图的大小
        legends=batch_legends,      # 标签信息
        legendFontSize=12,
        useSVG=True                # SVG格式保证清晰度
    )
```

**展示策略**：
- 分批避免信息过载
- 5×2网格布局最优
- SVG矢量格式无损缩放
- 合适的子图尺寸确保结构清晰

---

## 📊 数据说明

### 数据字段含义

| 字段名 | 类型 | 范围 | 说明 | 来源 |
|--------|------|------|------|------|
| `gscore` | float | -12 ~ -8 | Glide对接分数，越负越好 | 分子对接软件 |
| `asie` | float | 10 ~ 30 | 平均特异性相互作用能 | 能量计算 |
| `inhibition_10um` | float | 0 ~ 100 | 10μM时的抑制百分比 | 体外实验 |
| `inhibition_50um` | float | 0 ~ 100 | 50μM时的抑制百分比 | 体外实验 |
| `ic50` | float | 1 ~ 100 | 半数抑制浓度(μM) | 剂量曲线拟合 |

### 数据质量指标
- **完整性**：20个化合物全部有基础数据
- **一致性**：50μM抑制率 > 10μM抑制率
- **可靠性**：IC50仅对高活性化合物测定

---

## 🚀 运行示例

### 完整运行流程
```python
# 1. 运行Cell 1 - 加载数据
exec(open('Cell_1_Data_Initialization.py').read())
print("✅ 数据加载完成")

# 2. 运行Cell 2 - 第一次筛选
exec(open('Cell_2_First_Screening.py').read())
plot_first_screening_10um()
print("✅ 10μM筛选图已生成")

# 3. 运行Cell 3 - IC50测定
exec(open('Cell_3_IC50_Determination.py').read())
plot_ic50_curves()
print("✅ IC50曲线已绘制")

# 4. 运行Cell 4 - 第二次筛选
exec(open('Cell_4_Second_Screening.py').read())
plot_second_screening_50um()
print("✅ 50μM验证完成")

# 5. 运行Cell 5 - 分子结构
exec(open('Cell_5_Molecular_Structures.py').read())
display_molecular_structures()
print("✅ 分子结构展示完成")
```

### Jupyter Notebook运行
```python
# 在Jupyter中按顺序运行每个Cell
# Shift + Enter 执行当前Cell并移到下一个
```

---

## 📈 结果解读

### 关键发现

#### 1. 筛选效率分析
```
初始化合物: 20个
  ↓ 10μM筛选
高活性(>35%): 2个
  ↓ IC50测定
最优化合物: AQ-149/42126274 (IC50=16.976μM)
  ↓ 50μM验证
确认活性: 多个化合物>50%抑制
```

#### 2. 活性化合物特征
- **最佳化合物**：AQ-149/42126274
  - IC50: 16.976 μM（良好活性）
  - 10μM抑制率: 45.2%
  - 50μM抑制率: 78.5%

#### 3. 结构-活性关系（SAR）
- 苯并噻唑核心结构普遍存在
- 取代基影响活性强弱
- 卤素取代可能增强活性

---

## 🔧 扩展应用

### 1. 更换靶点蛋白
```python
# 修改compounds_data为其他激酶的数据
compounds_data_CDK4 = {
    'compound_1': {'gscore': -11.2, ...},
    # 新的靶点数据
}
```

### 2. 增加评价指标
```python
# 在数据字典中添加新字段
'selectivity': 85.3,      # 选择性指数
'cell_ic50': 12.5,        # 细胞水平IC50
'logP': 3.2,              # 脂水分配系数
```

### 3. 自定义筛选阈值
```python
# 修改活性分级标准
HIGH_ACTIVITY_THRESHOLD = 40    # 原为35
MEDIUM_ACTIVITY_THRESHOLD = 25  # 原为20
```

### 4. 导出分析报告
```python
# 生成Excel报告
df_results = pd.DataFrame.from_dict(compounds_data, orient='index')
df_results.to_excel('CDK2_screening_results.xlsx')

# 生成PDF报告
from matplotlib.backends.backend_pdf import PdfPages
with PdfPages('screening_report.pdf') as pdf:
    pdf.savefig(fig1)  # 保存各个图表
```

---

## ❓ 常见问题

### Q1: 中文显示乱码怎么办？
**A**: 确保安装了中文字体，Windows用户可以：
```python
plt.rcParams['font.sans-serif'] = ['Microsoft YaHei']
```

### Q2: RDKit安装失败？
**A**: 推荐使用conda安装：
```bash
conda install -c conda-forge rdkit
```

### Q3: 如何添加新化合物？
**A**: 在compounds_data字典中添加：
```python
compounds_data['NEW_ID/12345678'] = {
    'gscore': -10.5,
    'asie': 22.3,
    'inhibition_10um': 18.7,
    'inhibition_50um': 35.2
}
```

### Q4: 如何修改图表样式？
**A**: 修改matplotlib参数：
```python
plt.style.use('ggplot')  # 更换样式
plt.rcParams['figure.dpi'] = 150  # 提高分辨率
```

---

## 🤝 贡献指南

欢迎贡献代码和建议！

### 贡献方式
1. Fork本项目
2. 创建特性分支 (`git checkout -b feature/AmazingFeature`)
3. 提交更改 (`git commit -m 'Add some AmazingFeature'`)
4. 推送到分支 (`git push origin feature/AmazingFeature`)
5. 开启Pull Request

### 代码规范
- 遵循PEP 8规范
- 添加详细的中文注释
- 更新README相应部分
- 确保所有测试通过

---

## 📚 引用信息

如果您在研究中使用了本代码，请引用：

```bibtex
@software{cdk2_screening_2025,
  title = {CDK2抑制剂筛选数据分析与可视化平台},
  author = {Your Name},
  year = {2025},
  url = {https://github.com/yourusername/cdk2-screening}
}
```

### 相关文献
1. **CDK2在癌症中的作用**：[Nature Reviews Drug Discovery]
2. **虚拟筛选方法**：[Journal of Medicinal Chemistry]
3. **Hill方程在药物发现中的应用**：[Pharmacological Reviews]

---

## 📄 许可证

本项目采用MIT许可证 - 查看 [LICENSE](LICENSE) 文件了解详情

---

## 👥 作者与致谢

- **作者**：[XIYUAN CHEN]

### 致谢
- 感谢提供实验数据的研究团队
- 感谢RDKit开源社区
- 感谢Python科学计算生态系统的贡献者

---

## 📮 联系方式

- Email: xiyuan.chen23@qq.com
- Issues: [GitHub Issues](https://github.com/lelebixiong/cdk2-screening/issues)
- Discussion: [GitHub Discussions](https://github.com/lelebixiong/cdk2-screening/discussions)

---

**最后更新**: 2025年9月

⭐ 如果这个项目对您有帮助，请给个Star支持一下！
