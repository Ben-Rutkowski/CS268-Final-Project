Let $G = (\mathcal{N}, \mathcal{E})$ be a graph in which nodes $\mathcal{N}$ represents sites in which an EV charging sation can be build and edges $\mathcal{E}$ represents the connection of sites by drivable routes. Let $d:\mathcal{N} \times \mathcal{N} \rightarrow \mathbb{R}^{+}$ be the distance of road between the shorts route in $G$ between nodes $i$ and $j$. The parameters $f_i$, $F_i$ and $c_i$ represent the charging capacity, charging demand, and charging cost respectively, for each node $i$. Let $D$ be the average distance an EV can drive on a full battery without recharging. Define $\mathcal{N}_i^{\alpha D} = \{i\in \mathcal{N}|d(i,j)\leq \alpha D\}$ to be the set of nodes adjacent that are at most $\alpha D$ distance away from $i$, where $\alpha \in (0,1]$. Let $\hat{G} = (\hat{\mathcal{N}}, \hat{\mathcal{E}})$ where $\hat{\mathcal{N}} = \mathcal{N}$ and $\hat{\mathcal{E}} = \{(i,j)|i,j\in \mathcal{N}, d(i,j)\leq D, i\neq j\}.$

Let $0^i$ be a source node (not in $\mathcal{N}$) that provides flow directly to node $i$ and through the rest of the graph $\hat{G}$. Let $0\leq x_0^i \leq n$ (where $n = |\mathcal{N}|$) be the residue flow not consumed by the system. For each edge $(j,k)\in \hat{\mathcal{E}}$, let the flow on $(j,k)$ originating from source $0^i$ be represented by the variable $y_{jk}^i$. 

The design points of our objective function are boolean variables $x_i\in \{0,1\}$, which represent whether or not a charging station is build at node $i$, and continuous variables $y_{jk}^i$. The optimization problem is to minimize the cost function 
$
\begin{align}
    \text{minimize } \sum_{i=1}^n c_i x_i
\end{align}
$
with constraints
$
\begin{align}
    \sum_{j\in \mathcal{N}_i^{\alpha D}} f_j x_j & \geq F_i & \forall i \\
    x_i & = \{0, 1\} & \forall i \\
    x_0^i + y_{0i}^i & = n & \forall i \in \hat{\mathcal{N}} \\
    0 \leq y_{jk}^i & \leq n x_i x_k & \forall (j,k) \in \hat{\mathcal{E}} \cup (0^i, i), \forall i\in \hat{\mathcal{N}} \\
    \sum_{j|(j,k)\in \hat{\mathcal{E}}} y_{jk}^i & = x_i x_k + \sum_{l|(k,l)\in \hat{\mathcal{E}}} y_{kl}^i & \forall i,k \in \hat{\mathcal{N}} \\
    x_i \sum_{j\in \hat{\mathcal{N}}} x_j & = y_{0i}^i & \forall i\in \hat{\mathcal{N}} \\
    0 & \leq x_0^i & \forall i\in \hat{\mathcal{N}}.
\end{align}
$
