# Traffic Coordination Game

This repo contains Matlab simulation code for the Traffic Coordination Game. These simulations where developed as part of a group project in the course __Game Theory and Rationality ENM140__ at Chalmers University of Technology, Sweden, 2016.

This code may freely be used and built upon for private and educational purposes as long as the original code is not misrepresented. Proper disclosure would of course be appreciated.

## Game Description

(Due to lack of LaTeX support, or similar math formatting, this description might be less readable than it should be. Apologies.)

Consider a set of *N* players. All players wish to get from point A to point B, between which there are *m* paths. All players must then simultaneously choose a path, where the pay-off *P_i* for choosing path *i*, is given by

*P_i = m + 1 - i - c N_i*,

where *c* is a cost parameter and *N_i* is the number of __other__ players that chose the same path. The game is a single-round game with $m$ possible actions (the *m* paths) for every player.

-----
#### Author:
*Simon Nilsson (simnilss)*