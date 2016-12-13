PARTNO = undef;
zmin = -203.8e-3;
zmax = 196.2e-3;

R = 55e-3;

N = 128;
phis = [ for (i = [0:N-1]) i*360/N ];
avertices = [ for (phi = phis) [R*cos(phi), R*sin(phi), zmin] ];
bvertices = [ for (phi = phis) [R*cos(phi), R*sin(phi), zmax] ];

module radial_wall() {
    faces = concat([[N-1, 0, N, 2*N-1]],
        [ for (i = [0:N-2]) [i, i+1, i+N+1, i+N] ]);
    polyhedron(concat(avertices, bvertices), faces);
}

module acap() {
    verts = concat([[0, 0, zmin]], avertices);
    faces = concat([[0, N, 1]], [ for (i = [1:N-1]) [0, i, i+1]]);
    polyhedron(verts, faces);
}

module bcap() {
    verts = concat([[0, 0, zmax]], bvertices);
    faces = concat([[0, N, 1]], [ for (i = [1:N-1]) [0, i, i+1]]);
    polyhedron(verts, faces);
}

if (PARTNO == 0) {
    radial_wall();
} else if (PARTNO == 1) {
    acap();
} else if (PARTNO == 2) {
    bcap();
} else {
    radial_wall();
    acap();
    bcap();
}
