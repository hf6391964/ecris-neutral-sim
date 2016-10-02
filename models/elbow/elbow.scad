PARTNO = undef;
A = 4;//undef;
B = 3;//undef;

R = 1;

N = 128;
phis = [0 : 360/N : 360];
avertices = [ for (phi = phis) [A, R*cos(phi), R*sin(phi)] ];
bvertices = [ for (phi = phis) [R*cos(phi), B, R*sin(phi)] ];

module elbow() {    
    evertices = [ for (phi = phis) [R*cos(phi), R*cos(phi), R*sin(phi)] ];
    faces = [ for (i = [0:N]) [i, i+1, i+N+2, i+N+1] ];
    polyhedron(concat(avertices, evertices), faces);    
    polyhedron(concat(bvertices, evertices), faces);
}

module acap() {
    verts = concat([[A, 0, 0]], avertices);
    polyhedron(verts, [ for (i = [0:N]) [0, i, i+1] ]);
}

module bcap() {
    verts = concat([[0, B, 0]], bvertices);
    polyhedron(verts, [ for (i = [0:N]) [0, i, i+1] ]);
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
