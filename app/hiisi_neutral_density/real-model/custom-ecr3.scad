cylinder_length = 400;
cylinder_r_inner = 50;
cylinder_r_outer = 56;
groove_depth = 3;
groove_fillet = 3;
groove_width = 28;
injection_end_h = 4;
extraction_end_h = 4;

waveguide_lx = 20;
waveguide_ly = 10;

injection_tube_r = 6;
injection_tube_wall_thickness = 1;

injection_plug_h = 40;
injection_plug_r1 = 20;
injection_plug_r2 = 10;

grill_hole_r = 1.5;

$fa = 5;
$fs = 2;

module body() {
    difference() {
        cylinder(r=cylinder_r_outer,
            h=cylinder_length + injection_end_h + extraction_end_h);
        translate([0, 0, injection_end_h])
        linear_extrude(cylinder_length)
        offset(r=-groove_fillet) offset(delta=groove_fillet)
        offset(r=groove_fillet) offset(delta=-groove_fillet)
        union() {
            intersection() {
                union() {
                    circle(r=cylinder_r_inner);
                    for (i = [0:6]) {
                        rotate([0, 0, i/6 * 360])
                        resize([cylinder_r_inner + groove_depth, groove_width, 0])
                            translate([0, -0.5, 0]) square();
                    }
                }
                circle(r=cylinder_r_inner + groove_depth);
            }
        }
    }
}

module grill(r, R, N) {
    for (i = [0:N]) {
        rotate([0, 0, i/N*360])
            translate([R, 0, -1]) cylinder(r=r, h=injection_end_h+2);
    }
}

module assembly(preview) {
    difference() {
        body();
        if (!preview) {
            grill(grill_hole_r, 48, 75);
            grill(grill_hole_r, 44, 65);
            grill(grill_hole_r, 40, 60);
        }
    }
}

assembly(false);
