cylinder_length = 400;
cylinder_r_inner = 50;
cylinder_r_outer = 56;
groove_depth = 3;
groove_fillet = 3;
groove_width = 28;
injection_end_h = 4;
extraction_end_h = 4;

waveguide_lx = 10;
waveguide_ly = 20;
waveguide_wall_d = 2;
waveguide_h = 40;
waveguide_angle = 40;
waveguide_r = 25;
waveguide_fillet = 2;

injection_tube_r = 5;
injection_tube_wall_d = 1;
injection_tube_h = 40;
injection_tube_R = 47;

injection_plug_h = 40;
injection_plug_r1 = 20;
injection_plug_r2 = 10;

extraction_cyl_r = 35;
extraction_cyl_offset = 5;
extraction_cyl_h = 30;
extraction_hole_r = 4;

grill_hole_r = 1.5;
grill_sector_angle = 75;

$fa = 5;
$fs = 1;

cylinder_outer_length = cylinder_length + injection_end_h + extraction_end_h;

module body() {
    difference() {
        cylinder(r=cylinder_r_outer, h=cylinder_outer_length);
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
    for (k = [0:3]) {
        for (i = [0:N]) {
            rotate([0, 0, k*120 + (i/N - 0.5) * grill_sector_angle])
                translate([R, 0, -1]) cylinder(r=r, h=injection_end_h+2);
        }
    }
}

module extraction_cyl() {
    translate([0, 0, cylinder_outer_length - extraction_cyl_h]) difference() {
        cylinder(r=extraction_cyl_r, h=extraction_cyl_h);
        translate([0, 0, -1]) cylinder(r=extraction_hole_r, h=extraction_cyl_offset+2);
    }
}

module injection_tube() {
    rotate([0, 0, 60]) translate([injection_tube_R, 0, 0])
        cylinder(r=injection_tube_r+injection_tube_wall_d,
                 h=injection_tube_h);
}

module injection_tube_hole() {
    rotate([0, 0, 60]) translate([injection_tube_R, 0, -1])
        cylinder(r=injection_tube_r, h=injection_tube_h+2);
}

module injection_plug() {
}

module waveguide() {
    difference() {
        linear_extrude(waveguide_h) minkowski() {
            translate([waveguide_fillet + waveguide_r, 0, 0])
            resize([waveguide_lx + 2*waveguide_wall_d - 2*waveguide_fillet,
                waveguide_ly + 2*waveguide_wall_d - 2*waveguide_fillet, 0])
                translate([0, -0.5, 0]) square();
            circle(waveguide_fillet);
        }
        
        waveguide_dh = (waveguide_lx + 2*waveguide_wall_d) * tan(waveguide_angle);
        translate([waveguide_r, 0, waveguide_h - waveguide_dh])
            rotate([0, -waveguide_angle, 0])
            resize([2*(waveguide_lx + 2*waveguide_wall_d),
                    waveguide_ly + 2*waveguide_wall_d + 2, waveguide_dh + 2])
            translate([0, -0.5, 0]) cube();
    }
}

module waveguide_hole() {
    translate([waveguide_r + waveguide_wall_d, 0, -1])
        resize([waveguide_lx, waveguide_ly, waveguide_h + 2])
        translate([0, -0.5, 0]) cube();
}

module assembly(preview) {
    difference() {
        union() {
            body();
            extraction_cyl();
            for (i = [0:3]) {
                rotate([0, 0, i*120 + 60]) waveguide();
            }
            injection_tube();
        }
        
        for (i = [0:3]) {
            rotate([0, 0, i*120 + 60]) waveguide_hole();
        }
        
        injection_tube_hole();
        
        translate([0, 0,
            cylinder_outer_length - extraction_cyl_h + extraction_cyl_offset])
        cylinder(r=extraction_cyl_r-extraction_cyl_offset, h=extraction_cyl_h + 1);

        if (!preview) {
            grill(grill_hole_r, 46, 16);
            grill(grill_hole_r, 42, 14);
            grill(grill_hole_r, 38, 12);
        }
    }
}

assembly(false);
//assembly(true);
//projection(cut=true) translate([0, 0, -100]) assembly(true);
