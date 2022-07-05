# T-FSM: A Task-Based System for Massively Parallel Frequent Subgraph Pattern Mining from a Big Graph


## Serial Verision

Under the root directory of the project, execute the following commands to compile.
```
cd T-FSM/serial
make
```
Run serial FSM using the following command.
```
./run -file [GRAPHPATH] -freq [F] -maxNodes [MAXNODE]
```

GRAPHPATH: required, the input graph path

F: required, the user-give support threshold

MAXNODE: optional, parameter to set the maximum subgraph size (number of vertices)

Example:

Use the following command:
```
./run -file ../data/mico.lg -freq 9480
```

## Parallel Version
Under the root directory of the project, execute the following command to complie.
```
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

## Input

Input graph starts with 't N M' where N is the number of vertices and M is the number of edges. A vertex and an edge are formatted as 'v VertexID VertexLabel Degree' and 'e VertexID VertexID EdgeLabel' respectively. Note that we require that the vertex id starts from 0 and the range is [0,N - 1] where V is the vertex set. The following is an input sample
```
t 5 6
v 0 0 2
v 1 1 3
v 2 2 3
v 3 1 2
v 4 2 2
e 0 1 1
e 0 2 1
e 1 2 1
e 1 3 1
e 2 4 1
e 3 4 1
```
