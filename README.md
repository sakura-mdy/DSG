A7:
Our explanation is as follows: First, in each iteration, E_need represents the theoretical number of edges that need to be added in this iteration, and this value is independent of the specific way in which edges are added between nodes. Second, our edge addition process is carried out uniformly after the iteration is complete, following certain rules as described in A3. The final edge configuration is not unique, but the number of k-core followers, which we are concerned with, remains unchanged.

A8:
There was a typo in Algorithm 1. The iteration starts from the first one, and deg(u, k, 0) should be replaced with deg(u, k, 1).

A9:
It should indeed be v ∈ V(k, b, i-1) here.

A10:
We understand your concern regarding the conclusion drawn in Figure 2 based on experiments with only a single value of k. We would like to clarify that these experiments were intended as a motivating example to illustrate the initial behavior of our approach. The results are not meant to be definitive but to provide a preliminary understanding of the method's performance. In the appendix, we will provide further experiments with different parameter settings.

A11:
The examples in the paper did not further elaborate on the intuitive meaning of these concepts. I will revise this section to provide more clarity.

A12:
This was an oversight on our part. To maintain consistency, the λ-threshold component should be defined as a set of nodes with a dynamic degree less than or equal to k.

A13:
The way edges are inserted in line 17 follows the procedure described in A3. In line 6, the part mentioned in the pseudocode is converted into k-core followers, and the dynamic degrees of the nodes are updated accordingly.

A14:
As explained in my response in A3, no actual edges are inserted during the iteration process, and E_need is only a theoretical value. Therefore, it is possible to have two λ-shell components, each with an E_need of 0.5 (each having only one node with degree k-1). During the final edge addition stage, these two nodes can be connected to satisfy the conditions for both λ-shell components. This approach avoids the scenario in other methods where these two nodes are connected to nodes with degree greater than or equal to k.

A15:
Thank you for your suggestion. We will include additional experiments in the appendix to validate the results for varying k and b. It is worth noting that, theoretically, it has already been proven that our method outperforms previous approaches in terms of performance.
