PARTNO = undef;
A = 4;//undef;
B = 3;//undef;

R = 1;

N = 128;
phis = [ for (i = [0:N-1]) i*360/N ];
avertices = [ for (phi = phis) [A, R*cos(phi), R*sin(phi)] ];
bvertices = [ for (phi = phis) [R*cos(phi), B, R*sin(phi)] ];

module elbow() {
    evertices = [ for (phi = phis) [R*cos(phi), R*cos(phi), R*sin(phi)] ];
    faces = concat([[N-1, 0, 2*N, 3*N-1], [3*N-1, 2*N, N, 2*N-1]],
        [ for (i = [0:N-2]) [i, i+1, i+2*N+1, i+2*N] ],
        [ for (i = [0:N-2]) [i+2*N, i+2*N+1, i+N+1, i+N] ]);
    polyhedron(concat(avertices, bvertices, evertices), faces);    
}

module acap() {
    verts = concat([[A, 0, 0]], avertices);
    faces = concat([[0, N, 1]], [ for (i = [1:N-1]) [0, i, i+1]]);
    polyhedron(verts, faces);
}

module bcap() {
    verts = concat([[0, B, 0]], bvertices);
    faces = concat([[0, N, 1]], [ for (i = [1:N-1]) [0, i, i+1]]);
    polyhedron(verts, faces);
}

if (PARTNO == 0) {
    elbow();
} else if (PARTNO == 1) {
    acap();
} else if (PARTNO == 2) {
    bcap();
} else {
    elbow();
    acap();
    bcap();
}
