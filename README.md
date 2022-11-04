# T-FSM: A Task-Based System for Massively Parallel Frequent Subgraph Pattern Mining from a Big Graph


## Serial Verision

Under the root directory of the project, execute the following commands to compile.
```
cd T-FSM/serial
make
```
Run serial FSM using the following command. All frequent patterns are stored in grami_results.txt.
```
./run -file [GRAPHPATH] -freq [F] -maxNodes [MAXNODE]
```

GRAPHPATH: required, the input graph path

F: required, the user-given support threshold

MAXNODE: optional, parameter to constrain the maximum subgraph size (number of vertices)

Example:

Use the following command:
```
./run -file ../data/mico.lg -freq 9480
```

## Parallel Version
Under the root directory of the project, execute the following commands to complie. If the mined patterns are large (>=6 vertices), we highly suggest users enable #OPTIMIZED_MATCH option, which consumes slightly more memory. All frequent patterns are printed with #VERBOSE enabled.
```
cd T-FSM/parallel
make
```

Run parallel FSM using the following command.
```
./run -file [GRAPHPATH] -freq [F] -thread [T] -maxNodes [MAXNODE]
```

GRAPHPATH: required, the input graph path

F: required, the user-given support threshold

T: required, the number of threads, 32 by default

MAXNODE: optional, parameter to constrain the maximum subgraph size (number of vertices)

Example:

Use the following command:
```
./run -file ../data/mico.lg -freq 9480 -thread 4
```

## Fraction-Score Version
Under the root directory of the project, execute the following commands to complie.
```
cd T-FSM/fraction-score
make
```

Run Fraction-Score FSM using the following command.
```
./run -file [GRAPHPATH] -freq [F] -thread [T] -maxNodes [MAXNODE]
```

GRAPHPATH: required, the input graph path

F: required, the user-given support threshold, in terms of Fraction-Score metric

T: required, the number of threads, 32 by default

MAXNODE: optional, parameter to constrain the maximum subgraph size (number of vertices)

## Input

Input graph starts with 't N M' where N is the number of vertices and M is the number of edges. A vertex and an edge are formatted as 'v VertexID VertexLabel Degree' and 'e VertexID VertexID EdgeLabel' respectively. Note that we require that the vertexID starts from 0 and the range is [0,|V| - 1] where V is the vertex set. If the input graph is not edge-labelled, mark all edge labels as 1. The following is an input sample
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
