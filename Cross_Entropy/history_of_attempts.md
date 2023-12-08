# History of Attempts
Documenting the progression of research to include in the writeup

- Started with model that has boolean for each station
- Changed to model that has a single integer that represents the number of stations built
- Initially wanted to have a single distribution for whole vector
- But decided to have independent distributions for each station and a binomial dist for integer
- However, we want to eliminate ambiguity in each design point. We can represent the same arrangement of 
    stations by two different vectors (e.g. switching stations 1 and 2 in x). So we must correct each 
    vector after performing randomly generating one. We do this by sorting on the x-coord.

- We will find the top ten elite samples, force them to be feasible, but optimizing the cost, and then calculate the distrubution

- Seidel failed
- Implemented Simplex Method