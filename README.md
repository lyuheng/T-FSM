# T-FSM: A Task-Based System for Massively Parallel Frequent Subgraph Pattern Mining from a Big Graph


## Serial Verision

Execute the following commands to compile.
```
git clone https://github.com/lyuheng/T-FSM.git
cd T-FSM/serial
make
```
Run serial FSM using the following command.
```
./run -file [GRAPH_PATH] -freq [F] -maxNodes [MAXNODE]
```

GRAPH_PATH: required, the input graph path

F: required, the user-give support threshold

MAXNODE: optional, parameter to set the maximum subgraph size (number of vertices)

Example:

Use the following command:
```
./run -file ../data/mico.lg -freq 9480
```

## Parallel Version
Execute the following command to complie.
```
git clone https://github.com/lyuheng/T-FSM.git
cd T-FSM/parallel
make
```

Run serial FSM using the following command.
```
./run -file [GRAPHPATH] -freq [F] -thread [T] -maxNodes [MAXNODE]
```

GRAPHPATH: required, the input graph path

F: required, the user-give support threshold

T: required, the number of threads

MAXNODE: optional, parameter to set the maximum subgraph size (number of vertices)

Example:

Use the following command:
```
./run -file ../data/mico.lg -freq 9480 -thread 32
```
