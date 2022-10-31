# 有限元笔记

## 目录

[前言](#before)

[数学准备](#m_math)

&emsp;&emsp;[矩阵](#m_m_matrix)

&emsp;&emsp;[张量](#m_tensor)

&emsp;&emsp;[线性空间](#m_linear)

[离散系统下的数学模型](#m_stand)

&emsp;&emsp;[静态问题](#jump_to_stand)

[连续系统下的数学模型](#m_continue)

&emsp;&emsp;[变分法](#jump_to_variational)

&emsp;&emsp;[变分法下的近似解](#jump_to_variational_usolve)

----

## <span id='before'>前言</span>
本文档旨在通过学习有限元原理，以及使用程序（目前先决定用Python，后面可能会使用matlab与Fortran，方便我熟悉它们的编写，此外这三种语言不能说是分门别类，只能说是痛苦面具，面向对象、面向过程（我都不确定Fortran算不算严格意义上的面向过程，或者干脆面向块得了）和脚本语言的激烈碰撞……）
相关代码将会转移到我的github项目中，后续会添加github对应的地址，本文中所有代码块使用的均为python语言。

 
----


## <span id="m_math">数学准备</span>

### <span id='m_m_matrix'>矩阵</span>（随时更新）
在这个部分，将会说明一些矩阵上要注意的事，尽管如此，一些最基本的矩阵知识将不记录，像是*特征值*，*正定矩阵*，*稀疏矩阵*，*对角矩阵*等基础的矩阵定义则不再阐述，这里主要收录一些不怎么常见，或是在编程中需要用到的矩阵。

**对称矩阵**

对于一个矩阵$A_{n \times n}$来说，若有$A_{ij}=A_{ji},i=1,2,...,n,j=1,2,...,n$，则称该矩阵为对称矩阵。

该矩阵也可由下式定义：
$$
A=A^T
$$
对于这样的一个矩阵，若完整地把它完全存储起来，随着矩阵维数的扩张，吞噬的内存将相当竞人，这也会增加计算成本，因此可以使用**上三角矩阵**来代替一个对称矩阵并将其储存起来，而对于大多数有限元问题，**稀疏矩阵**是被常常用到的，这意味着矩阵中大多数的元素都为0。

这样即便用一个上三角矩阵来储存似乎也没能省下多少空间（事实上如何实现只存一个上三角就是一个好问题），因此这里使用一个**一维数组**来储存上面的信息。

在此之前，先介绍另外一个矩阵：*带状矩阵*。对于一个带状矩阵$A$来说，其元素$a_{ij}$满足：

$$
a_{ij}=0, \text{when }j>i+m_A
$$

其中$m_A$表示半带宽，$2m_A+1$则是带宽，如下例所示为一个带宽为5，半带宽为2的对称带状矩阵：

$$
A=
\left[\begin{array}{l}
3 & 2 & 1 & 0 & 0 \\
2 & 3 & 2 & 1 & 0 \\
1 & 2 & 3 & 2 & 1 \\
0 & 1 & 2 & 3 & 2 \\
0 & 0 & 1 & 2 & 3
\end{array}\right]
$$
对于这个对称带状矩阵来说，可以确定的是在带宽外的元素会始终为0，而带宽内的元素随着计算的进行可能为0也可能不为0。

因此我们可以通过一种方法来构造一个一维数组$l$来储存上述的对称带状矩阵$A$：

$$
l(1)=a_{11},l(2)=a_{22},l(3)=a_{21},l(4)=a_{33},l(5)=a_{32},l(6)=a_{31},...
$$

下面是使用代码描述该方法的过程

	# m为半带宽, n为阶数, A为已经定义好的一个nxn对称带状矩阵
	m = 2
	n = 5
	l = []
	for i in range(n):
		l.append(A[i][i])
		for j in range(i - 1, i - m - 1, -1):
			if j < 0:
				continue
			l.append(A[i][i - 1])

虽然看起来这样不容易做计算，但实际上十分方便。

### <span id='m_tensor'>张量</span>（随时更新）
对于空间中一向量$u$来说，对于其基向量$e_i$，有：
$$
u=\sum{u_ie_i}
$$
其中$u_i$为各基上的坐标，使用张量可简单表示为：

$$
u=u_ie_i
$$

可见右下角的指标就算换成$j$等表示完全不影响最后的结果，因此称这样的指标为**哑指标**。

克罗内克符号：

$$
\delta_{ij}=
\begin{cases}
\delta_{ij}=1, i=j \\
\delta_{ij}=0, i \neq j
\end{cases}
$$

###  <span id="m_linear">线性空间</span>

在这部分将会简单地对**线性空间**这个概念进行介绍。但先让我们把这个抽象概念变的具体一点。

在高中的时候我们都学过集合，而线性空间即是**满足一些特定规则**的向量集合（因此有时也被成为向量空间），但这里的向量并不是那个熟悉的，有方向有大小的那个向量$\overrightarrow{a}$，这个向量更加的广义，或者说，在线性空间中的所有元素（不管它是一个标量还是矢量还是函数）都被称为向量。

关于它的判断方法，我并不想过多赘述——因为这不是我们现在要关心的，反而重要的是由它引申出的两个概念：线性相关和线性无关。

> 毕竟要判断一个向量集合是不是一个线性空间，你要一条一条地检查共计8条规则。

在学习矩阵时，其实我们已经接触到线性相关和线性无关这两个概念了，在线性空间中，这个概念几乎可以看成当时我们学习这个概念时的延展，这里就不抛砖引玉了，直接说结论：

若对于线性空间$V$，有向量组$\{k_i|k_i\in V,i=1,2,...,n\}$判断他们是否线性相关的依据即为下式是否存在非0解：

$$
\sum_{i=1}^nk_ix_i=0
$$

若不存在非0解，则说明向量组$\{k_1,k_2,...,k_n\}$是线性无关的，同理反之则是线性相关。

补充一点，若向量组$\{k_1,k_2,...,k_n\}$线性无关，且$V$中任一元素都可由$\{k_1,k_2,...,k_n\}$线性表示：

$$
\sum_{i=1}^nk_ix_i=A,\quad A\in V
$$

则我们将向量组$\{k_1,k_2,...,k_n\}$称为是$V$的一组基，而上式的解$\{x_1,x_2,...,x_n\}$称为向量$A$在以$\{k_1,k_2,...,k_n\}$为基下的坐标。

现在让我们把这个问题变得抽象一点，要多抽象呢？

假设现在有一个线性空间$M$，其中的所有向量都是一个线性空间$N_i$，当然这些向量之间存在一些关系，这些我们并不关心，其中$\{N_1,N_5,N_9\}$是线性无关的，那么我们可以得到个什么结论呢？类比上面的式子：

$$
N_1x_1+N_5x_2+N_9x_3=0
$$

这个式子只存在0解，可能会有疑问：这里作为向量的线性空间和一个标量的乘积，线性空间之间的加法该如何定义？不用担心这个，总会有一个定义的，但通过上式可以说明：**这三个向量无法通过线性变换彼此表示**。

这个概念非常重要，因为在后文会见到这些概念。

----

### <span id='m_stand'>离散系统下的数学模型</span>

在很多工程问题中，需要先建立数学模型，再使用诸如有限元等计算方法进行求解。而有限元作为一种离散的方法，需要对不同的数学模型进行分类讨论。

这里的离散系统指的是*集中参数模型*，而连续介质力学模型则在后面的连续系统的数学模型中提到。

在一个离散系统中，为了描述系统的状态，往往需要如下几步：

 1. 离散化问题，即把一个问题理想化为单元（们）的组合体
 2. 平衡单元：确定单元平衡要求
 3. 单元组装：通过单元之间的连接来建立位置状态量的联立方程组
 4. 计算响应：通过求解联立方程得到状态变量，并由单元见的连接要求，计算每个单元的响应。

#### <span id='jump_to_stand'>静态问题</span>

所谓静态问题也就是系统相应不随时间变化而变化。

更简单来说就是在《弹性力学》中遇到的种种问题，或是《材料力学》中的一些问题，在一个杆件系统或者平面问题上施加一个静态荷载，整个系统同时也是静态的。更一般的，对于系统中的每个状态分量$U_i$都有：

$$
\frac{\partial U_i}{\partial t}=0
$$
上式说明各个分量与时间均无关，同时这种问题也是最简单的一类问题之一，因此作为有限元的切入点还是蛮不错的。

接下来用一个很简单的例子来说明这一过程：

<center>
<a href="https://sm.ms/image/yg6mw8ZtGELTcKY" target="_blank">
<img src="https://s2.loli.net/2022/10/25/yg6mw8ZtGELTcKY.jpg" />
<br>
简单的线性单元（过于简陋了）
</a>
</center>

上面的图是个很简陋的线性单元结构，两端节点分别有$U_1$和$U_2$的位移相应以及$F_1$和$F_2$的载荷，显然单元$k_1$的平衡方程可写为：

$$
\left[
\begin{array}{}
k_1 & -k_1 \\
-k_1 & k_1
\end{array}
\right]
\left[
\begin{array}{l}
U_1 \\
U_2
\end{array}
\right]=
\left[
\begin{array}{l}
F_1 \\
F_2
\end{array}
\right]
$$

通过上面的矩阵式可以很看出矩阵第一行的系数决定了节点1对于单元$k_1$的贡献，同理第二行则是节点2对单元$k_1$的贡献。

在此基础上对于更为复杂的线性离散单元系统进行整合：

<center>
<a href="https://sm.ms/image/WtaMpmiqQPnb49I" target="_blank">
<img src="https://s2.loli.net/2022/10/25/WtaMpmiqQPnb49I.jpg" >
<br>
<center>
摘自《有限元法理论、格式与求解方法》
</center>
</a>
</center>

在上面的系统中，有3个节点和5个单元，这里仅以单元2的平衡方程为例：

$$
\left[
\begin{array}{}
k_2 & -k_2 & 0 \\
-k_2 & k_2 & 0 \\
0 & 0 & 0
\end{array}
\right]
\left[
\begin{array}{}
U_1 \\
U_2 \\
U_3
\end{array}
\right]=
\left[
\begin{array}{l}
F_1^{(2)} \\
F_2^{(2)} \\
F_3^{(2)}
\end{array}
\right]
$$

而整个系统的平衡方程经整理后可写作：

$$
\left[
\begin{array}{}
k_1+k_2+k_3+k_4 & -(k_2+k_3) & -k_4 \\
-(k_2+k_3) & k_2+k_3+k_5 & -k_5 \\
-k_4 & -k_5 & k_4+k_5
\end{array}
\right]
\left[
\begin{array}{l}
U_1 \\
U_2 \\
U_3
\end{array}
\right]=
\left[
\begin{array}{l}
R_1 \\
R_2 \\
R_3
\end{array}
\right]
$$

上式中的系数矩阵我们可以记作$K$，也即**刚度矩阵**。同时上式也可写为：

$$
KU=R
$$

对于刚度矩阵有：

$$
K=\sum_{i=1}^5{K^{(i)}}
$$

显然，刚度矩阵为一个对称矩阵，而且在多数实际情况下，它是一个稀疏的**对称带状矩阵**，这样就可以用我们在[矩阵](#m_m_matrix)那儿提到的应对对称带状矩阵的存储方式，在实际的编程中进行存储。

这样直接了当写出整个刚度矩阵并直接求解的方法也被成为**直接刚度法**，同时在线性系统中（如刚才的例子），所有刚度矩阵中的元素都是**常数(const)**，当然这不意味着在非线性系统中无法使用直接刚度法，只不过现在这部分内容还不涉及到我要整理的内容，因此**暂且略过**。（谁知道以后会不会补上呢）

同时别忘了这是静态系统，也就是当受力平衡后，**整个系统的状态变量组成的泛函应在鞍点**，即对由$U_i(i=1,2,...,n)$组成的泛函$\Pi (U_1, ... U_n))$应有：

$$
\delta \Pi = 0
$$

而

$$
\delta \Pi = \frac{\partial \Pi}{\partial U_1}\delta U_1+\frac{\partial \Pi}{\partial U_2}\delta U_2+...+\frac{\partial \Pi}{\partial U_n}\delta U_n=0
$$

因此一定有：

$$
\frac{\partial \Pi}{\partial U_i}=0
$$

因为在上式中，对于状态变量的变分不一定为0（除非在边界上），因此为保证泛函变分为0才有上式的结论。

而对于诸如一维杆的问题、二维平面问题等经典的力学问题（包括结构力学的一些桁架问题等），该泛函可表示为：

$$
\Pi = \boldsymbol{u}-\boldsymbol{w}
$$

这里$\boldsymbol{u}$表示系统的应变能，$\boldsymbol{w}$表示载荷的总势能。

对于符合胡克定律的线性系统，上式又可表达为：

$$
\boldsymbol{u}=\frac{1}{2}ku^2, \boldsymbol{w}=Pu
$$

对泛函取变分，这样我们能得到这样的结果：

$$
\delta \Pi = \frac{\partial \Pi}{\partial u}\partial u = 
(ku-P)\partial u,\frac{\delta \Pi}{\partial u}=ku-P
$$

由前文知：

$$
\frac{\delta \Pi}{\partial u}=0\Rightarrow ku=P
$$

可以发现，使用变分法得到的平衡方程与直接刚度法得到的平衡方程是一致的。使用矩阵来表示这个过程：

$$
\boldsymbol{u}=\frac{1}{2}U^TKU, \boldsymbol{w}=U^TR
$$


将其代入泛函表达式，注意到$KU=R$有：

$$
\Pi=-\frac{1}{2}U^TR
$$

此即为系统平衡时的位移值，也即为所求状态变量的值，因为此时泛函$\Pi$有极小值——这也是我们要取得的目的之一。使用变分法优点是能够针对问题快速建立平衡方程，缺点就是其物理意义解释起来比直接刚度法要困难一些。

## <span id='m_continue'>连续系统下的数学模型</span>

在刚刚我们介绍完离散系统的数学模型，通过变分法或者直接刚度法我们可以直观地获得系统的平衡方程，但对于连续系统中我们发现直接刚度法变得没有那么好用了——因为所有“点”都是连续的，这样该如何解决平衡方程呢？

这里同样是两种不同的方法，分为微分法和变分法。

对于微分法不过多阐述，因为跟我目前所要解决的问题关系并不大，这里主要说一下变分法。

> 这里仅对微分法做最基本的介绍，简单来说微分法是通过一个微分式：
> $A(x,y)\frac{\partial^2u}{\partial x^2}+2B(x,y)\frac{\partial^2y}{\partial x \partial y}+C(x,y)\frac{\partial^2 u}{\partial y^2}=\varphi(x,y,u,\frac{\partial u}{\partial x},\frac{\partial u}{\partial y})$
> 来确定系统的平衡方程，其中通过对ABC的不同情况，可分为
> $B^2-AC=\begin{cases}<0, 椭圆 \\ =0,抛物线 \\ >0,双曲线\end{cases}$
> 不同的边界条件来解决实际问题，这类方程主要用于特征值的求解问题中，目前不在我的研究领域中，因此这里简单跳过。

### <span id='jump_to_variational'>变分法</span>

这里依旧使用上一节中的泛函$\Pi$，设泛函中的状态变量的最高阶导数为m，控制微分方程的最高导数为2m阶*，称此类问题为$C^{m-1}$变分问题，对于$C^{m-1}$变分问题的本质边界，大多数产生于对应的指定的位移和转角，本质边界条件中，$C^{m-1}$变分问题的导数最多为(m-1)阶。

> *笔者注：我觉得还是按泛函中状态量的最高阶导数理解比较好理解，举个最简单的例子，在一维杆的问题中，很明显状态变量$u$在沿杆方向只有一阶导是有意义的（应变），单纯的二阶段不存在物理意义（欧拉梁理论中根据平衡方程$\frac{\partial}{\partial x}(EA\frac{\partial u}{\partial x})=P$也可以看出在平衡方程中是存在二阶导的），所以该问题是一种$C^0$变分问题。

而对应自然边界，也即力边界问题，在结构力学*中自然边界条件对应于指定的边界力和边界力矩，所以这些边界的最高阶导数为m至2m-1。
> *笔者注：这里可能是因为翻译问题，应该是一般的力学边界都称之为自然边界，或是力边界。

在继续往下之前，需要知道有这几种运算方法，假设每个给定的$x$的函数$F$依赖于$v$（状态变量），$dv/dx$, ..., $d^pv/dx^p$，其中$p=1,2,...$。函数$F$一阶变分定义如下：

$$
\delta F=\frac{\partial F}{\partial v}\delta v+\frac{\partial F}{\partial(dv/dx)}\delta\frac{dv}{dx}+...
$$

之后我们构造一个函数$\varepsilon\eta(x)$，其中$\varepsilon$为常数，$\eta(x)$为任一光滑的函数，其值在对应的本质边界条件上为0，同时有$\eta(x)=\delta v(x)$。对其求$n$阶导就会变成：

$$
\frac{d^n\eta}{dx^n}=\frac{d^n\delta v}{dx^n}=\delta(\frac{d^nv}{dx^n})
$$

这样对于$F$的变分定义式又可以写为：

$$
\delta F=\lim_{\varepsilon\rightarrow0}\frac{F[v+\varepsilon\eta,\frac{d(v+\varepsilon\eta)}{dx},...,\frac{d^p(v+\varepsilon\eta)}{dx^p}]-F(v,\frac{dv}{dx},...,\frac{d^pv}{dx^p})}{\varepsilon}
$$

其和导数的定义形式类似，同时某些运算形式也和导数的运算性质高度相似：

$$
\begin{array}{}
\delta(F+Q)=\delta F+\delta Q \\
\delta(FQ)=(\delta F)Q+F(\delta Q) \\
\delta(F)^n=nF^{n-1}\delta F
\end{array}
$$

### <span id='jump_to_variational_usolve'>变分法的近似解</span>

上一节已经提出了对于一个连续系统中状态变量的泛函变分形式，当然对于一些已知泛函变分形式的问题，可能可以通过分部积分等方法得到精确的泛函解（这要求不低的数学素养），但更多时候，我们并不需要也不能得到一个问题的精确解，这时候就需要通过种种方法得到一个近似解来尽可能地得到我们想要的结果。

我们采用如下的微分形式来描述一个系统稳态问题的状态：

$$
L_{2m}[\varphi]=r
$$

其中，$L_{2m}$为线性微分算子，$\varphi$是状态变量，$r$是强迫函数，该问题的解应满足边界函数：

$$
B_i[\varphi]=q_i|_{S_i}, i=1,2,3,...
$$

> 笔者注：这里之所以取$L_{2m}$为微分算子的符号，是和之前的$C^m$变分问题保持一致，例如当该问题为$C^0$变分问题时，在平衡微分方程中最高阶导数为2阶，此时$m=1$。

其中对于微分算子$L_{2m}$，特别要注意当它为正定算子时，满足对称条件（关于为什么我们要注意它是否正定，后面会给出讨论）：

$$
\int_D(L_{2m}[u])vdD=\int_D(L_{2m}[v])udD
$$

以及正定性：

$$
\int_D(L_{2m}[u])udD>0
$$

下面用一个简单的一维欧拉梁来举例子。

对于一个一维杆的稳态响应，采取如下微分方程计算：

$$
-EA\frac{\partial^2u}{\partial x^2}=0
$$

满足边界条件：

$$
u|_{x=0}=0;\quad EA\frac{\partial u}{\partial x}|_{x=L}=R
$$

通过比较对应的算子，我们可以得到如下结论：

$$
L_{2m}=-EA\frac{\partial^2}{\partial x^2};\quad \varphi=u;\quad r=0
$$

同样的，对于边界条件通过对比可得：

$$
\begin{array}{ll}
B_1=1; & q_1=0 \\
B_2=EA\frac{\partial }{\partial x}; & q_2=R
\end{array}
$$

接下来我们确定算子$L_{2m}$是否正定，这时我们可以抛去一切外载荷的影响，只考虑结构本身。因此我们取$R=0$：

$$
\begin{array}{lll}
\displaystyle\int^L_0{-EA\frac{\partial^2u}{\partial x^2}vdx} & = & \displaystyle-EA\frac{\partial u}{\partial x}v|^L_0+\int_0^LEA\frac{\partial u}{\partial x}\frac{\partial v}{\partial x}dx \\
 & = & \displaystyle-EA\frac{\partial u}{\partial x}v|^L_0 + EAu\frac{\partial u}{\partial x}|^L_0-\int_0^LEA\frac{\partial^2v}{\partial x^2}udx
\end{array}
$$

别忘了我们的假设，这时上式中的$\displaystyle\frac{\partial u}{\partial x}$在两个边界上均为0，因此上式经过简化后写为：

$$
\displaystyle\int^L_0{-EA\frac{\partial^2u}{\partial x^2}vdx}=\int_0^L-EA\frac{\partial^2v}{\partial x^2}udx
$$

满足了正定算子的对称性，其实完全可以确定它也是正定的，因为有：

$$
\int^L_0-EA\frac{\partial^2u}{\partial x}udx=\int^L_0EA\left(\frac{\partial u}{\partial x}\right)^2dx>0
$$

接下来，让我们来看看两种方法：**加权余量法和Ritz法**来求解我们需要的状态变量，这两种办法都会把解的基本形式假设为：

$$
\overline{\varphi}=\sum^n_{i=1}a_if_i
$$

其中，$f_i$是彼此线性无关*的试函数，$a_i$是待定系数。

> *笔者注：详见[数学准备-线性空间](#m_linear)

这里考虑加权余量法，其中余量定义式：

$$
R=r-L_{2m}\left[\sum_{i=1}^na_if_i\right]
$$

对于精确解，该余量当然等于0，但明显这种情况很难出现，因此我们尽力要做的是让余量尽量逼近于0，在这里我们争取让余量的加权平均为0。这里我们介绍两种方法：伽辽金（Galerkin）法和最小二乘法：

**Galerkin法**：在该方法中，参数$a_i$由$n$个方程确定：

$$
\int_DRdD=0;\quad i=1,2,...,n
$$

**最小二乘法**：该方法中，余量平方的积分关于$a_i$是最小的：

$$
\frac{\partial}{\partial a_i}\int_DR^2dD=0;\quad i=1,2,...,n
$$

经过变换后：

$$
\int_DRL_{2m}[f_i]dD=0;\quad i=1,2,...,n
$$

还有两种方法分别是配点法和子域法，这两种方法由于实现起来可能会相当困难，因此一般只考虑Galerkin法和最小二乘法。

接下来介绍里茨（Ritz）法，相比于最小余量法，Ritz法的核心是将上面提到的式子$L_{2m}[\varphi]=r$使用一个等价的泛函$\Pi$表示，同时对于试函数$\overline\varphi$也代入$\Pi$中，并让$\Pi$的驻值为0：$\delta\Pi=0$，得出关于参数$a_i$的$n$个联立式子：

$$
\frac{\partial\Pi}{\partial a_i}=0;\quad i=1,2,...,n
$$

在这其中，试函数$f_i$的选择相对放松，只需要满足本质边界条件而无需去验证是否满足自然边界条件（因为自然边界条件包含在了泛函$\Pi$中）当对应的微分算子$L_{2m}$是正定的时候，泛函$\Pi$的极小值即为最小值。

> 当然，找到一个满足自然边界条件的试函数能更精确地解答问题——但这往往是很困难的，因此一般都会选择一个满足本质边界条件的试函数来得到一个近似解。

接下来举一个例子来说明Ritz法：

现有一悬臂杆，共长180cm，左侧100cm部分较细，横截面积为1$cm^2$，右侧80cm部分较粗，横截面积随长度变化而变化，变化函数为：$S=(1+\frac{x-100}{40})^2cm^2$，坐标原点为最左端固定端，杆长沿$x$正向，右端自由端受集中力$R=100N$，方向沿$x$正向，该结构的总势能如下：

$$
\Pi=\int_0^{180}\frac{1}{2}EA\left(\frac{\partial u}{\partial x}\right)^2dx-100u|_{x=180}
$$

现在进行试函数的假设，这里分两种情况讨论：

（1）假设试函数为

$$
u=a_1x+a_2x^2
$$

同时让泛函$\Pi$的变分为0，即：

$$
\delta\Pi=\int_0^{180}EA\left(\frac{\partial u}{\partial x}\right)\delta\left(\frac{\partial u}{\partial x}\right)dx-100\delta u|_{x=180}
$$

该式经过分部积分后得到：

$$
\begin{array}{c}
\displaystyle\frac{d}{dx}\left(EA\frac{du}{dx}\right)=0 \\
\\
\displaystyle EA\frac{du}{dx}|_{x=180}=100
\end{array}
$$

同时还有个边界条件：$u|_{x=0}=0$，利用这些边界条件可知：

$$
\begin{array}{c}
\displaystyle u=\frac{100}{E}x,\quad 0 \leqslant x \leqslant 100 \\
\\
\displaystyle u=\frac{10000}{E}+\frac{4000}{E}-\frac{4000}{E\left(1+\frac{x-100}{40}\right)},\quad 100 \leqslant x \leqslant 180
\end{array}
$$

同时，杆中的应力精确解应为：

$$
\begin{array}{c}
\displaystyle \sigma=100;\quad 0 \leqslant x \leqslant 100 \\ \\
\\
\displaystyle \sigma=\frac{100}{1+\frac{x-100}{40}};\quad 100 \leqslant x \leqslant 180
\end{array}
$$

通过一个简单的程序，我们可以把上图画出来（假设有$E=10^6$，代码可以在[这里]()找到）： 

 <center>
<a href="https://sm.ms/image/5RFtHDizgG1l7rb" target="_blank"><img src="https://s2.loli.net/2022/10/31/5RFtHDizgG1l7rb.png"  width="300" high="300"></a>
</center>

以及对应的应变图：

<center>
<a href="https://sm.ms/image/MYGyD5pehB9R3OC" target="_blank"><img src="https://s2.loli.net/2022/10/31/MYGyD5pehB9R3OC.png" width="300" high="300"></a>
</center>

对于这些
