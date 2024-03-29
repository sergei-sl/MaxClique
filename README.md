# Description

Realisation of branch and bound algorithm for finding exact solution for the Maximum clique problem using greedy clique heuristic for lower bound on each step.<br>
Project also contains Cliquer solver:  https://users.aalto.fi/~pat/cliquer.html (can be used for debugging purposes)
This solver can use Domain Filtering (currently disabled):  https://link.springer.com/chapter/10.1007/3-540-45749-6_44
This is C++ realization for branch-and-bound algorithm: https://www.researchgate.net/publication/222747193_A_Branch-and-Bound_Algorithm_for_the_Maximum_Clique_Problem
In every branch we solve the corresponding LP relaxation for the maximum clique problem using Cplex
## Instructions

To run the executable you should write in the command line:
- [executable_name] [path_to_the_file] [time_limit_in_seconds]

Default timeout is 0

## Test results (for those graphs, where clique was found before the program was terminated by 3600 sec timeout)

Graph name|#Nodes|#Edges|Found clique size|Time (ms)
---|---|---|---|---
brock200_2.clq|200|9876|12|3145298
brock200_1.clq|200|14834|21|2867141
brock200_4.clq|200|13089|17|2048040
c-fat200-1.clq|200|1534|12|25143
c-fat200-2.clq|200|3235|24|461
c-fat200-5.clq|200|8473|58|24496
c-fat500-2.clq|500|9139|26|23966
c-fat500-5.clq|500|23191|64|4
hamming6-2.clq|64|1824|32|354
hamming6-4.clq|64|704|4|0
hamming8-4.clq|256|20864|16|25911
johnson8-2-4.clq|28|210|4|0
johnson8-4-4.clq|70|1855|14|70
johnson16-2-4.clq|120|5460|8|4780
MANN_a9.clq|45|918|16|912
p_hat300_1.clq|300|10933|8|70
p_hat300_2.clq|300|21928|25|637892
p_hat500_1.clq|500|31569|9|325
san200_0.7_1.clq|200|13930|30|2231
san200_0.7_2.clq|200|13930|18|60428
san200_0.9_1.clq|200|17910|70|170
san200_0.9_2.clq|200|17910|60|163
san200_0.9_3.clq|200|17910|44|19648
san400_0.5_1.clq|400|39900|13|1206612
C125.9.clq|125|6963|34|226118
gen200_p0.9_44.clq.txt|200|17910|44|500
gen200_p0.9_55.clq.txt|200|17910|55|1290

