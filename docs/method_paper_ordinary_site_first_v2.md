# 方法：结构约束下的表面组合突变筛选

本文方法用于在给定蛋白结构上选择少量、可解释且计算稳定性风险受控的表面双点或三点突变体，为后续功能实验提供扰动探针。方法将位点组合的几何支持、侧链化学扰动、表位先验和稳定性预测分开计算，不把任一代理指标直接解释为抗体结合亲和力或实验表位。筛选流程依次包括结构输入与残基映射、表面区域划分、单点替换准入、位点集合枚举、侧链机制排序、三原则队列评估、SPURS 组合筛选、FoldX 重复确认、正式资格复核和精简面板选择。

筛选参数在候选评估前固定。表面暴露分类、预测结构置信度分档和 FoldX 稳定性变化幅度采用相关文献或官方方法说明中的判据作为参考；其在本流程中的准入方向和组合判定规则另行明确。几何尺度、SPURS 分档、表位百分位和覆盖增益等采用本研究预设值。各项参数的来源与适用范围列于补充表 S1。

## 1. 结构输入、残基映射与质量控制

输入包括分析序列及其坐标范围、选定的 PDB/mmCIF 目标链、结构类型（实验结构或预测结构）、链编号和可用的结构置信度。对预测结构，保留模型文件与对应置信度来源；AlphaFold 3 生成的模型按其预测置信度解释，不视为实验结构[4]。序列位置、结构 author 编号及插入码、规范化预测链位置同时保存。所有突变定义均使用“序列位置—野生型残基—目标残基”的明确对应关系；野生型不一致、位置重复或超出序列范围的替换不进入候选空间。

表面暴露以残基总溶剂可及面积（solvent-accessible surface area，SASA）相对于该残基最大允许面积的比例表示：

\[
\mathrm{rSASA}_i=\frac{\mathrm{SASA}_i}{\mathrm{MaxASA}(a_i)}.
\]

本文以 rSASA≥0.25 作为表面分析准入阈值。25% 的界限参考已有溶剂可及性研究中划分暴露与埋藏残基的做法[14]；本文以 Tien 等的最大允许面积归一化[2]，因此沿用的是分类界限，而非原研究的全部归一化设置。SASA 采用 Shrake–Rupley 数值积分和 1.4 Å 溶剂探针计算[1]，实现使用 Biopython 的 Shrake–Rupley 接口[3]。ProCaliper 作为另一实现用于 SASA/rSASA 一致性核对[5]，不构成独立方法验证，也不作为有限尺寸抗体探针可达性的替代计算。

对预测结构，残基 pLDDT<70 或缺少相应置信度时标记为结构复核项。70 的界限参考 AlphaFold 数据库的置信度解释：70–90 通常支持较可靠的主链预测，50–70 属于低置信度区间[15,17]。本文据此设置结构复核门；pLDDT≥70 不单独保证侧链构象或复合物界面准确。实验结构的 B-factor 不当作 pLDDT 使用。结构复核、表面暴露和突变许可分别记录，任一坐标可用并不意味着该位点可以进入正式推荐。

## 2. 表面区域划分

在蛋白表面残基图上构建连通区域。每个显示区域须满足以下几何条件：区域内部连通；最大表面测地直径为 30 Å；最大欧氏直径为 22 Å；区域面积不少于 50 Å²；至少包含 4 个残基。区域规模按本次分析中可映射的残基总数 L 设置软目标：

| 分析残基数 L | 目标区域大小 | 允许的偏好范围 |
| --- | --- | --- |
| L<120 | 7 个残基 | 6–8 |
| 120≤L≤150 | 8 个残基 | 6–10 |
| L>150 | 12 个残基 | 10–15 |

目标区域大小和区域数量是设计偏好，连通性、直径、面积和最少残基数是硬性几何检查。区域提供候选归属和区域覆盖率分母；全局精简面板不替代逐区域推荐。

## 3. 单点替换空间与序列保护

对每个已映射且通过表面准入的位点，枚举除野生型外的 19 种标准氨基酸替换。替换进入组合位点池须同时满足：位点位于结构图、已被标记为可分析表面、无残基级硬排除、无结构复核项，并具有有限且有效的 SPURS 单点稳定性预测。SPURS 是基于蛋白生成模型表征的稳定性预测方法[8]，其输入表征与 ESM2、ProteinMPNN 等蛋白语言模型和序列设计模型有关[9,10]。

单点预测值 s 按单向上限分层，数值方向统一为去稳定化 ΔΔG（kcal/mol）：

| SPURS 单点值 s | 分类 | 是否进入组合候选池 |
| --- | --- | --- |
| s≤0.5 | PREFERRED | 是，仍须满足其他条件 |
| 0.5<s≤1.0 | ACCEPTABLE | 是，仍须满足其他条件 |
| 1.0<s≤2.0 | REVIEW | 否 |
| s>2.0 | REJECT | 否 |
| 缺失或非有限值 | NOT_SCORED | 否 |

判断使用 s 本身而非 |s|。0.5、1.0 和 2.0 的 SPURS 分档为本研究预设的风险筛选标准；SPURS 原文[8]提供预测方法依据，不作为这些分界值已获验证的依据。单点分数只用于构造候选空间，组合突变仍须单独评分。

序列保护规则包括：保护野生型 C、G、P，禁止自动引入 C、G、P；保留二硫键、已注释的糖基化位点和翻译后修饰；保留其他残基级硬排除。对每个多点组合，在完整参考序列上比较替换前后的 N-X-S/T 糖基化序列子集合；X 不能为 P。组合只有在不新增、不删除既有序列子且位置、野生型残基和重复位点均合法时才进入后续计算。N-X-S/T 序列子规则依据天冬酰胺连接糖基化的共识序列[13]。

## 4. 位点集合枚举与几何支持

对每个区域 R，收集至少存在一个合格单点替换且位于结构图上的残基，枚举其中全部二元和三元位点子集。双点集合的两个位点须在表面图上可达，且最短测地距离不超过 24 Å；三点集合要求三对位点均满足同一条件。几何量只依赖位点集合，不依赖具体氨基酸替换，因此同一位点集合的所有替换共享几何结果。

设 a_j 为表面图中残基 j 的合格映射面积，d_G 为表面图最短路径距离。对半径 r∈{9,12,15} Å，定义位点集合 S 的邻域并集：

\[
U_r(S)=\bigcup_{i\in S}\{j:d_G(i,j)\leq r\}.
\]

区域覆盖率及三尺度平均覆盖率为：

\[
C_r(S\mid R)=
\frac{\sum_{j\in U_r(S)\cap R}a_j}
{\sum_{j\in R}a_j},
\qquad
\overline C(S\mid R)=
\frac{C_9+C_{12}+C_{15}}{3}.
\]

对每一尺度，边际面积定义为该集合相对于最佳删去一个位点后所增加的邻域面积：

\[
M_r(S)=
\max\left[
0,\ A(U_r(S))-\max_{i\in S}A(U_r(S\setminus\{i\}))
\right],
\]

其中 A(T)=Σ_{j∈T}a_j。三个尺度的边际面积取平均作为组合边际支持。邻域重叠度先计算各位点邻域之间的面积 Jaccard 比值，再对位点对和三个尺度分别取平均。直接位点面积为 S 中各位点合格映射面积之和。

三点集合还须通过区域增益门。设 D* 为同一区域内覆盖率最高的几何合法双点集合，则：

\[
g_R(T)=\overline C(T\mid R)-\overline C(D^*\mid R),
\qquad g_R(T)\geq 0.05.
\]

基线 D* 不要求是 T 的子集，也不要求已通过稳定性评分；缺少三点覆盖、区域面积或双点基线时，增益视为未知，三点集合不能越过该门。该规则用于防止三点组合仅因位点数增加而在区域覆盖上缺乏实质收益。

## 5. 侧链机制描述与候选排序

每个实际替换由六个侧链化学轴描述：名义电荷、体积、疏水性、氢键供体数、氢键受体数和芳香性。对含 k 个替换的集合 S，轴 x 的聚合变化为：

\[
D_x(S)=\frac{1}{k}\sum_{i=1}^{k}
\left|x(a_i^{\mathrm{mut}})-x(a_i^{\mathrm{WT}})\right|.
\]

暴露轴使用各位点侧链暴露比例的最小值：

\[
e_i=
\frac{\mathrm{SASA}^{\mathrm{context}}_{\mathrm{sidechain},i}}
{\mathrm{SASA}^{\mathrm{isolated\ residue}}_{\mathrm{sidechain},i}},
\qquad
e_{\min}(S)=\min_{i\in S}e_i.
\]

侧链 SASA 采用 Shrake–Rupley 算法、1.4 Å 探针、每原子 960 个采样点和固定的元素半径，只计算重原子。侧链不完整、残基类型无对应描述符或孤立残基分母缺失时保留为未知；任一位点缺少某轴时，该组合的相应轴也记为未知，不以零填补。含组氨酸的替换，其名义电荷、氢键供体和受体变化保留为未知。电荷和疏水性另报告带符号净变化及方向抵消比例 \(O_x=1-|\sum_i\Delta x_i|/\sum_i|\Delta x_i|\)，但排序使用绝对变化均值。

机制比较同时考虑暴露最小值与六个绝对变化轴。只有缺失轴模式相同的两个候选才进行支配判断：若 A 在所有可比轴上不低于 B，且至少一个轴更高，则 A 支配 B。排序先按证据完整性和缺失模式分组，再在组内按 Pareto 层比较。

为形成可解释的机制优先级，六个变化轴和暴露轴被归入三组：暴露、电荷/极性、体积/芳香。每个轴在本次调用参考分布的有限唯一值中进行序数归一化；若某轴有 n>1 个不同值，则第 i 个值的效用为 i/(n−1)，从 0 开始编号。常数轴不提供区分效用。每组取该组可用轴效用的最大值，再将三组效用从小到大排序形成 leximin 向量，优先最大化最低组效用，随后依次比较其余分量。参考分布按当前调用范围建立，不能把归一化效用当作跨蛋白的绝对尺度。

候选的最终排序依次考虑证据完整性、缺失模式、Pareto 层、leximin 向量，再以几何支持、风险和稳定候选标识解决平局。该排序用于在有限候选空间内确定验证顺序，不等同于对全部可能替换组合的全局最优搜索。

## 6. 每位点代表替换与位点集合预排

每个位点按三种角色选择真实替换代表：电荷/极性、体积/芳香和综合机制。角色内排序依次比较证据完整性、角色效用、较低的单点 SPURS 值和氨基酸字母序。初始代表优先覆盖三种角色；若候选耗尽且尚未获得合格推荐，则按同一排序逐层补充，单个位点最多保留 7 种代表。分层扩展用于控制组合空间，所有结构资格、单点准入、三点增益、SPURS 和 FoldX 门均保持不变。

为避免在全部替换组合上立即展开，每个位点集合先用少量实际替换组成见证组合。见证组合须通过完整序列检查；机制队列用同一机制排序函数比较见证组合，并以其中最优见证项预排位点集合。见证组合只是队列排序代理，不是已完成的 SPURS/FoldX 验证结果，也不保证代表该集合全部替换的最优解。

## 7. 三条独立设计原则及排序

每个区域分别建立表位支持、机制扰动和区域覆盖三条独立队列。每条队列在首个满足正式科学门和展示条件的候选处完成；同一构建体被不同原则选中时合并原则标签，但不重复计为多个构建体。

### 7.1 表位支持原则

表位支持原则适用于抗原分析。每个位点使用已评分的表面表位百分位；局部峰或百分位不低于 0.80 的位点标记为受支持。局部峰定义为该位点分数高于表面图上的全部直接邻居；上游若提供局部峰标记则沿用该标记。表位先验可由结构型 B 细胞表位倾向模型产生，本文采用 DiscoTope-3.0 作为软先验[6]。0.80 和局部峰规则为本研究的相对优先级设置，不对应 DiscoTope-3.0 原文给出的固定性能阈值；原文讨论了分数校准及抗原长度相关偏差[6]。该分数表示表位倾向，不是特定抗体的结合概率；构象表位预测方法的性能边界已有系统评述[7]。

设 \(p_i\) 为位点 i 的表位百分位，集合 S 的表位排序键为：

\[
E_1(S)=\frac{\sum_{i\in S}\mathbf{1}(i\ \text{受支持})}{|S|},
\qquad
E_2(S)=\frac{\sum_{i\in S}p_i}{|S|}.
\]

候选先按 E1 降序，再按 E2 降序。集合中任一位点缺少有限表位先验时，该集合在表位队列中记为证据缺失并跳过；其他原则仍可独立运行。Binder 分析不使用抗原表位先验，该原则记为不适用。

### 7.2 机制扰动原则

位点集合先按见证组合的机制排序预排；实际访问集合后，对物化替换使用第 5 节的暴露、Pareto 和 leximin 排序。该原则偏好具有明确侧链化学变化且暴露支持充分的组合，不以单一稳定性分数替代机制解释。所有候选仍须通过相同的 SPURS 和 FoldX 阈值。

### 7.3 区域覆盖原则

位点集合依次按以下顺序排序：区域平均覆盖率更高；平均边际面积更大；平均邻域重叠更低；直接位点面积更大。具体替换出现平局时，再比较初始风险、表位分数、突变阶数和稳定候选标识。该原则优先选择空间互补、能在区域内形成较完整表面扰动的组合。

### 7.4 双点优先与三点提前访问

一般情况下先访问双点集合，再访问三点集合。三点集合只有在通过区域增益硬门，并且在当前原则下优于其精确双点亲本时，才允许提前访问。精确双点亲本由见证三点中删除一个替换得到，其余位点和目标氨基酸保持不变。

| 原则 | 三点提前访问的比较条件 |
| --- | --- |
| 机制扰动 | 见证三点按机制排序支配最强的精确双点亲本 |
| 区域覆盖 | 三点平均区域覆盖至少比最强精确双点亲本高 0.05 |
| 表位支持 | 表位证据完整，且见证三点的表位排序键优于可比的精确双点亲本 |

未提前访问但通过区域增益门的三点集合仍可在双点集合之后进入验证。队列按提前访问标记、突变阶数、原则优先级、区域增益和位点集合标识依次排序；只对当前队首作完成判定。

## 8. SPURS 组合筛选与 FoldX 重复确认

组合评分使用实际完整的双点或三点突变签名，不把单点分数相加作为组合结果。SPURS 组合值 \(s_{\mathrm{multi}}\) 的分层为：

| SPURS 组合值 | 分类 | 是否进入 FoldX |
| --- | --- | --- |
| \(s_{\mathrm{multi}}\leq 1.0\) | PASS | 是 |
| \(1.0<s_{\mathrm{multi}}\leq 2.0\) | REVIEW | 否 |
| \(s_{\mathrm{multi}}>2.0\) | OVER_LIMIT | 否 |

双点和三点采用相同阈值。只有 PASS 候选进入 FoldX。FoldX 使用经验力场评估突变对蛋白稳定性的影响[11,12]，并与 SPURS 形成方法来源不同的复核证据。

FoldX 的去稳定化上限设为 2.0 kcal/mol。该幅度参考 Juritz 等使用 FoldX 预测时以 ±2 kcal/mol 区分近中性与稳定性扰动替换的先例[16]。原研究采用双向分类且分析单点替换；本文采用单向去稳定化上限，保留低于负向界限的预测稳定化组合，并将其用于双点和三点候选的设计筛选。因此，2.0 kcal/mol 是有文献先例的筛选幅度，不代表组合突变保持实验功能的通用界限。

为降低单次采样波动对判定的影响，FoldX 采用本研究预定义的重复采样规则；该规则并非文献[16]的采样方案。设首轮三个有效值为 \(f_0,f_1,f_2\)，并令

\[
h=\sum_{t=0}^{2}\mathbf{1}(f_t>2.0).
\]

| 首轮结果 | 后续判定 |
| --- | --- |
| 三值均有效且 h=0 | 直接 PASS |
| 三值均有效且 h≥2 | SCIENTIFIC_FAILURE |
| 三值均有效且 h=1 | 扩展至共 10 个有效槽，按十槽中位数判定 |
| 槽缺失、无效或身份不合格 | 记录证据异常，不以部分中位数放行 |

十槽判定要求完整、有限且编号对应的 10 个配对 ΔΔG 值，并通过

\[
\operatorname{median}(f_0,\ldots,f_9)\leq 2.0\ \mathrm{kcal/mol}
\]

才记为确认通过。十个值的中位数取排序后第 5、6 个值的平均；阈值 2.0 本身属于通过范围。首轮只有一个值超过 2.0 时，即使前三槽中位数不超过 2.0，也必须完成十槽确认。通过十槽判定的候选保留采样敏感性标记，不因中位数通过而删除实际跨门事实。

## 9. 正式科学门与区域推荐

候选通过 FoldX 后，按 candidate_id 独立复核正式资格，至少包括：不触及保护规则；位点和突变组合合法且连通；没有结构复核项；原生结构与规范化结构的表面暴露分类未发生改变；区域连通性、距离和直径检查有效；单点证据完整且无硬排除；SPURS 与 FoldX 均有有限、完整且通过的组合证据；三点集合满足最终区域覆盖增益。模型评分可按完整突变签名复用，但区域上下文、结构复核和资格状态不能跨 candidate_id 共享。

每个区域的每条适用原则在首个满足正式科学门和展示条件的候选处结束。不同原则选中同一组合时合并原则标签；不同组合分别保留。逐区域推荐不受全局精简面板预算裁剪。未完成队列须记录明确原因，包括证据缺失、没有合法候选、候选耗尽、技术失败或达到运行预算；这些状态不能与“已获得合格推荐”混同。

## 10. 全局精简面板与可辨识性

在逐区域推荐之外，另行选择有限构建体面板。贪心选择首先比较新增直接突变位点的合格表面积，其次比较几何边际支持、邻域重叠、假设信息增益和风险。默认最多选择 2 个高于双点阶数的主面板构建体；总面板规模由具体任务预算控制。相同完整构建体不重复选择；达到最小构建体数后，如没有正的直接覆盖增益，则停止。

对假设 \(h_a,h_b\) 和已选面板中的构建体 c，保留实际假设响应 \(r(c,h)\)，并计算：

\[
d(h_a,h_b)=
\sqrt{
\frac{\sum_{c\in\mathrm{panel}}
\left[r(c,h_a)-r(c,h_b)\right]^2}{8}
}.
\]

分母固定为 8，不随实际面板构建体数改变。\(d<0.1\) 的假设对记为近似不可区分。该距离用于报告面板设计中的未覆盖假设和难区分分组，不构成实验可分辨性的证明。

## 11. 解释边界

本方法输出的是在给定结构、表面几何和计算稳定性模型下可优先实验验证的突变组合。rSASA、测地距离、表面覆盖、侧链暴露和机制描述符均为结构代理量；SPURS 和 FoldX 为预测稳定性证据，不是实测表达量、折叠状态、细胞表面展示或抗体结合亲和力。表位倾向分数只用于软排序，不代表特定抗体的结合概率。

候选空间由有限替换代表、几何距离、序列保护和稳定性阈值共同限制，因此方法不保证遍历全部可能组合，也不保证找到全局最优解。文献支持的分类界限与本研究自设参数的适用范围见补充表 S1；若用于不同蛋白、不同结构置信度或不同实验体系，应在外部分析中报告阈值敏感性。历史任务仍按其保存的协议解释，当前方法定义不追溯改写既往结果。

## 补充表 S1. 筛选阈值及其依据

| 参数 | 本文取值或规则 | 依据及适用范围 |
| --- | --- | --- |
| 表面准入 | rSASA≥0.25 | 参考 25% 暴露/埋藏分类界限[14]；归一化参考值采用 Tien 等[2]；等于阈值时纳入为本文约定。 |
| 预测结构复核 | pLDDT<70 或置信度缺失 | 70 参考 AlphaFold 置信度分档[15,17]；缺失时复核及据此限制正式推荐为本文规则。 |
| 区域几何与规模 | 测地直径≤30 Å、欧氏直径≤22 Å、面积≥50 Å²、残基数≥4；按 L 采用 7/8/12 的软目标及第 2 节范围 | 本研究预设，用于限制区域空间跨度和规模；无直接支持这些具体数值的文献依据。 |
| 位点集合与邻域 | 位点对测地距离≤24 Å；邻域半径 9、12、15 Å | 本研究预设的局部组合与多尺度几何支持参数；不解释为实测抗体接触半径。 |
| 三点增益 | 区域覆盖增益≥0.05；覆盖原则提前访问时相对精确双点亲本亦须增加≥0.05 | 本研究预设，均指覆盖率增加至少 5 个百分点；两个比较的双点基线不同，见第 4、7 节。 |
| 单点与组合 SPURS | 单点分档 0.5/1.0/2.0；单点准入≤1.0；组合准入≤1.0，1.0–2.0 为复核区间 | 本研究预设；文献[8]支持模型来源，未直接验证本流程的分档和准入阈值。 |
| 表位支持 | 表面百分位≥0.80 或局部峰 | 本研究预设的相对排序规则；模型与校准依据见[6]，不赋予该规则原模型固定阈值的灵敏度或特异度。 |
| FoldX 稳定性 | 去稳定化 ΔΔG≤2.0 kcal/mol | 2 kcal/mol 幅度有 FoldX 单点研究先例[16]；本文采用单向上限，并用于组合设计。重复采样与中位数确认规则由本研究设定。 |
| 假设距离 | 固定分母 8；d<0.1 为近似不可区分 | 本研究预设的设计诊断尺度，未经实验可辨识性校准。 |

每位点代表替换上限和主面板三点构建体上限属于本研究的搜索与面板规模约束，不具有文献阈值的含义。上述文献先例也不构成本方法参数已在目标抗原或实验体系中完成校准的证据。

## 参考文献

1. Shrake A, Rupley JA. Environment and exposure to solvent of protein atoms. Applications to lysozyme and insulin. *J Mol Biol*. 1973;79(2):351–371. <https://doi.org/10.1016/0022-2836(73)90011-9>
2. Tien MZ, Meyer AG, Sydykova DK, Spielman SJ, Wilke CO. Maximum allowed solvent accessibilities of residues in proteins. *PLoS One*. 2013;8(11):e80635. <https://doi.org/10.1371/journal.pone.0080635>
3. Cock PJA, Antao T, Chang JT, et al. Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics*. 2009;25(11):1422–1423. <https://doi.org/10.1093/bioinformatics/btp163>
4. Abramson J, Adler J, Dunger J, et al. Accurate structure prediction of biomolecular interactions with AlphaFold 3. *Nature*. 2024;630(8016):493–500. <https://doi.org/10.1038/s41586-024-07487-w>
5. Rozum JC, Ufford H, Im AK, et al. ProCaliper: functional and structural analysis, visualization, and annotation of proteins. *Bioinformatics Advances*. 2025;5(1):vbaf275. <https://doi.org/10.1093/bioadv/vbaf275>
6. Høie MH, Gade FS, Johansen JS, et al. DiscoTope-3.0: improved B-cell epitope prediction using inverse folding latent representations. *Front Immunol*. 2024;15:1322712. <https://doi.org/10.3389/fimmu.2024.1322712>
7. Cia G, Pucci F, Rooman M. Critical review of conformational B-cell epitope prediction methods. *Brief Bioinform*. 2023;24(2):bbac567. <https://doi.org/10.1093/bib/bbac567>
8. Li Z, Luo Y. Generalizable and scalable protein stability prediction with rewired protein generative models. *Nat Commun*. 2026;17:891. Published online 20 December 2025. <https://doi.org/10.1038/s41467-025-67609-4>
9. Lin Z, Akin H, Rao R, et al. Evolutionary-scale prediction of atomic-level protein structure with a language model. *Science*. 2023;379(6637):1123–1130. <https://doi.org/10.1126/science.ade2574>
10. Dauparas J, Anishchenko I, Bennett N, et al. Robust deep learning–based protein sequence design using ProteinMPNN. *Science*. 2022;378(6615):49–56. <https://doi.org/10.1126/science.add2187>
11. Guerois R, Nielsen JE, Serrano L. Predicting changes in the stability of proteins and protein complexes: a study of more than 1000 mutations. *J Mol Biol*. 2002;320(2):369–387. <https://doi.org/10.1016/S0022-2836(02)00442-4>
12. Schymkowitz J, Borg J, Stricher F, et al. The FoldX web server: an online force field. *Nucleic Acids Res*. 2005;33(Web Server issue):W382–W388. <https://doi.org/10.1093/nar/gki387>
13. Kornfeld R, Kornfeld S. Assembly of asparagine-linked oligosaccharides. *Annu Rev Biochem*. 1985;54:631–664. <https://doi.org/10.1146/annurev.bi.54.070185.003215>
14. Petersen B, Petersen TN, Andersen P, Nielsen M, Lundegaard C. A generic method for assignment of reliability scores applied to solvent accessibility predictions. *BMC Struct Biol*. 2009;9:51. <https://doi.org/10.1186/1472-6807-9-51>
15. Varadi M, Bertoni D, Magana P, et al. AlphaFold Protein Structure Database in 2024: providing structure coverage for over 214 million protein sequences. *Nucleic Acids Res*. 2024;52(D1):D368–D375. <https://doi.org/10.1093/nar/gkad1011>
16. Juritz E, Fornasari MS, Martelli PL, Fariselli P, Casadio R, Parisi G. On the effect of protein conformation diversity in discriminating among neutral and disease related single amino acid substitutions. *BMC Genomics*. 2012;13(Suppl 4):S5. <https://doi.org/10.1186/1471-2164-13-S4-S5>
17. EMBL-EBI, Google DeepMind. AlphaFold Protein Structure Database: Frequently asked questions, “How confident should I be in a prediction?”. <https://www.alphafold.ebi.ac.uk/faq>. Accessed 4 November 2026.
