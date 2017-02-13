PARTNO = undef;

injection_R = 48;
injection_h = -17.3;
injection_angle = 30;
injection_r = 4.5;
injection_N = 50;

surround_cyl_dz = -450;
surround_cyl_r = 75;
surround_cyl_h = 510;

module chamber() {
    import("stl-external/plasmannakemat2-trimmed.stl", convexity=20);
}

module surrounding_cylinder() {
    translate([0, 0, surround_cyl_dz])
        cylinder(r=surround_cyl_r, h=surround_cyl_h);
}

module circle_surface(r, n) {
    verts = concat([[0, 0, 0]],
        [ for (i = [0:n]) [r*cos(i*360/n), r*sin(i*360/n), 0]]);
    faces = concat([[0, n, 1]], [ for (i = [1:n-1]) [0, i, i+1]]);
    polyhedron(verts, faces);
}

module injection_surface() {
    translate([injection_R*cos(injection_angle),
        injection_R*sin(injection_angle), injection_h])
        circle_surface(injection_r, injection_N);
}

module demo_injection_surface() {
    translate([injection_R*cos(injection_angle),
        injection_R*sin(injection_angle), injection_h - 1])
        circle(injection_r, $fn=injection_N);
}

translate([0, 0, -203.8]) mirror([0, 0, 1]) {
    if (PARTNO == 1) {
        chamber();
    } else if (PARTNO == 2) {
        surrounding_cylinder();
    } else if (PARTNO == 3) {
        injection_surface();
    } else {
        color([0, 0, 1, 0.1]) chamber();
        color([1, 0, 0, 0.5]) injection_surface();
        color([0, 1, 0, 0.05]) surrounding_cylinder();
    }
}
