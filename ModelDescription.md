# Electric Vehical Charging Station Location Optimization

### Basic Overview
The variables in our model are the locations and capcity of electric vehical charging stations under the constraint that all of the stations must meet the charging demand of every site in the area. Here sites may refer to neighborhoods, cities, counties, etc. depending on the scale of the application. The demand of a site may be the amount of EV's in the area. The capacity of a station may be the amount of chargers present in the station. We want to minimize the cost building all of the stations. 

### Model Description
The fixed parameters of this problem are the demand sites. We will call them *nodes* in this paper. There are $N$ total fixed nodes in the model. Each node has coordinates $(u_i, v_i)$ and a charging demand $F_i$ for all $i=1,\dots,N$. The other fixed parameters of this model are the initial cost of building a station $c$ and the cost of scaling the capacity of a station $f.$

The variables in this problem are properities of each of the *stations*. The design variable $\mathbf{x}$ is stored as the vector
$$
\mathbf{x} = 
\begin{pmatrix}
    x_1 \\ y_1 \\ z_1 \\ b_1 \\ \vdots \\ x_M \\ y_M \\ z_M \\ b_M
\end{pmatrix}
$$
where $x_k, y_k$ are the coordinates of the position $k$-th station and $z_k$ is the capacity of the $k$-th station. The variable $b_k \in \{0, 1\}$ is a boolean value that represents whether or not station $k$ is built.

### Cost Function
We want to minimize the cost of building all of the stations. The cost to build station $k$ is given by the formula $(c + fz_k)b_k.$ Here we scale the capacity cost $f$ by the capacity of the station $z_k$ and add $c$ the intial cost of building. We multiply the entire function by $b_k$ to indicate if the station is built (If the station is not built, $b_k = 0$ and the cost of the station is $0$). The maximum number of stations is $M$. Thus our design variable $\mathbf{x}\in \mathbb{R}^{4M}.$ 

The total cost to build all stations it thus given by 
$$
\begin{align}
    C(\mathbf{x}) & =\sum_{k=1}^{M} (c+fz_k)b_k
\end{align}
$$

### Constraints
The stations must have position and capacity that meet the demands of every node. A node has its demand met if all of the stations have capacity great enough to meet it. However, the further away a station is to a node, the less effective it's capacity is to meeting the demand of a node. This follows an inverse square relation. That is, if node $i$ and station $k$ are distance $d_{ik}$ apart, the supply given to node $i$ from station $k$ is given by $z_k / {d_{ik}}^2.$ Note that ${d_{ik}}^2 = (u_i - x_k)^2 + (v_i - y_k)^2.$ Thus, the total supply of that node $i$ recieves from every station is then given by
$$
\begin{align}
    S_i(\mathbf{x}) = \sum_{k=1}^M \frac{z_k}{(u_i - x_k)^2 + (v_i - y_k)^2}b_k
\end{align}
$$
(we multiply by $b_k$ per station once again because only built stations can provide supply).

In order for $\mathbf{x}$ to be feasible we want each node's demand $F_i$ to be met. The supply of node $i$ must be greater than or equal to it's demand i.e. $S_i(\mathbf{x}) \geq F_i.$

Lastly, we require that $b_k \in \{0, 1\}.$ In order to implement this formally on a continous $\mathbf{x}$, we require that ${b_k}^2 - b_k = 0.$ The left hand side is a polynomial whose only roots are $0$ and $1$. Thus $b_k \in \{0, 1\} \Leftrightarrow {b_k}^2 - b_k = 0.$ However, when optimizing in practice, we will not get exact values of $0$ and $1$. This is not a problem because each formula is multiplied by a $b_k$ and the difference between muliplying by $0$ and almost $0$ or by $1$ and almost $1$ is negligable in the limit.

### Formal Problem
We have a maxium of $M$ stations to be built yeilding a $\mathbb{R}^{4M}$ optimization problem. We want to determine the location and capacity of up to $M$ stations that minimizes the total cost of building the stations while meeting the demands of every node. Formally we represent this as the following.
$$
\begin{align}
    \text{minimize} \quad   & C(\mathbf{x}) = \sum_{k=1}^{M} (c+fz_k)b_k \\ 
    \text{subject to} \quad & g_i(\mathbf{x}) =  F_i - \sum_{k=1}^M \frac{z_k}{(u_i - x_k)^2 + (v_i - y_k)^2}b_k \leq 0 & \forall i\in 1,\dots,N \\
    & h_k(\mathbf{x}) = {b_k}^2 - b_k = 0 & \forall k = 1,\dots, M
\end{align}
$$