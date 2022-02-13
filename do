#!/bin/bash

################
### Compile code
################
g++ -std=c++0x -O2 ip3.cc -o ip3


######################
### Square wave 1 cell
######################
./ip3 -dt .025 -DT 20 -d .0 .0 -s 1 -T 600 -ca 0 0 -w 0.0 -prob 0.0  -kbath 4  -ca_sqONOFF 1 -ca_sqAMP .004 -ca_sqP 20000 data/square/1.sp -o  >data/square/1.hst 2>data/square/1.ca&
./ip3 -dt .025 -DT 20 -d .0 .0 -s 1 -T 300 -ca 0 0 -w 0.0 -prob 0.0  -kbath 4  -ca_sqONOFF 1 -ca_sqAMP .006 -ca_sqP 6666  data/square/2.sp -o  >data/square/2.hst 2>data/square/2.ca&
./ip3 -dt .025 -DT 20 -d .0 .0 -s 1 -T 300 -ca 0 0 -w 0.0 -prob 0.0  -kbath 4  -ca_sqONOFF 1 -ca_sqAMP .010 -ca_sqP 4000  data/square/3.sp -o  >data/square/3.hst 2>data/square/3.ca&
./ip3 -dt .025 -DT 20 -d .0 .0 -s 1 -T 300 -ca 0 0 -w 0.0 -prob 0.0  -kbath 4  -ca_sqONOFF 1 -ca_sqAMP .014 -ca_sqP 2857  data/square/4.sp -o  >data/square/4.hst 2>data/square/4.ca


############
### 2 cells
############
./ip3 -dt .025 -DT 1 -d .3 .3 -s 2 -iapp 1.25 1.25 -T 1000 -w 0.006 -prob 1 -kbath 8 -gleak 3.35 3.35 -fsca 0.02 0.02 data/2cell/1.sp -o  >data/2cell/1.hst 2>data/2cell/1.ca&
./ip3 -dt .025 -DT 1 -d .3 .3 -s 2 -iapp 2.00 2.00 -T 400  -w 0.006 -prob 1 -kbath 8 -gleak 3.35 3.35 -fsca 0.03 0.03 data/2cell/2.sp -o  >data/2cell/2.hst 2>data/2cell/2.ca&
./ip3 -dt .025 -DT 1 -d .3 .3 -s 2 -iapp 2.75 2.75 -T 400  -w 0.006 -prob 1 -kbath 8 -gleak 3.35 3.35 -fsca 0.06 0.06 data/2cell/3.sp -o  >data/2cell/3.hst 2>data/2cell/3.ca&
./ip3 -dt .025 -DT 1 -d .3 .3 -s 2 -iapp 3.50 3.50 -T 400  -w 0.006 -prob 1 -kbath 8 -gleak 3.35 3.35 -fsca 0.10 0.10 data/2cell/4.sp -o  >data/2cell/4.hst 2>data/2cell/4.ca

###########
## Network
###########
./ip3 -dt .025 -DT 20 -d .3 .3 -s 400 -fc 1.0 1.0 -T 300 -w 0.15 -prob 0.13 -kbath 5.5  -fsca 0.06 0.06 	data/net/1.sp -o  >data/net/1.hst 2>data/net/1.ca&
./ip3 -dt .025 -DT 20 -d .3 .3 -s 400 -fc 1.0 1.0 -T 300 -w 0.15 -prob 0.13 -kbath 6.5  -fsca 0.07 0.07 	data/net/2.sp -o  >data/net/2.hst 2>data/net/2.ca&
./ip3 -dt .025 -DT 20 -d .3 .3 -s 400 -fc 1.0 1.0 -T 300 -w 0.15 -prob 0.13 -kbath 7.5  -fsca 0.085 0.085 	data/net/3.sp -o  >data/net/3.hst 2>data/net/3.ca&
./ip3 -dt .025 -DT 20 -d .3 .3 -s 400 -fc 1.0 1.0 -T 300 -w 0.15 -prob 0.13 -kbath 8.5  -fsca 0.0975 0.0975   	data/net/4.sp -o  >data/net/4.hst 2>data/net/4.ca







