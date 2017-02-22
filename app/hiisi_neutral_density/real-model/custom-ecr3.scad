PARTNO = undef;

injection_end_h = 6;
extraction_end_h = 4;

waveguide_lx = 7.9;
waveguide_ly = 15.8;
waveguide_wall_d = 1;
waveguide_h = 54 + injection_end_h;
waveguide_angle = 60;
waveguide_r = 32;
waveguide_fillet = 1;

injection_tube_r = 5.25;
injection_tube_wall_d = 1;
injection_tube_h = 60 + injection_end_h;
injection_tube_R = 48;

injection_plug_h = 31;
injection_plug_r1 = 32;
injection_plug_r2 = 14.6;
injection_plug_button_r = 10;
injection_plug_button_h = 6;
injection_plug_plate_d = 6;
injection_plug_plate_h = 33;
injection_plug_plate_r = 42;
injection_plug_plate_clearance = 3;
injection_plug_bottom_r = 37.3;
injection_plug_bottom_h = 15;
injection_plug_h_total = injection_plug_plate_h;

extraction_cyl_r = 39;
extraction_cyl_offset = 1.9;
extraction_cyl_h = 30;
extraction_hole_r = 4;

cylinder_length = 397 + injection_plug_h_total + extraction_cyl_h;
cylinder_r_inner = 50;
cylinder_r_outer = 56;
groove_depth = 4.5;
groove_fillet = 5;
groove_width = 20;

grill_hole_r = 1.5;
grill_sector_angle = 80;

$fa = 10;
$fs = 1.5;

cylinder_outer_length = cylinder_length + injection_end_h + extraction_end_h;

module body() {
    difference() {
        cylinder(r=cylinder_r_outer, h=cylinder_outer_length);

        translate([0, 0, injection_end_h])
        linear_extrude(cylinder_length)
        offset(r=-groove_fillet) offset(delta=groove_fillet)
        offset(r=groove_fillet) offset(delta=-groove_fillet)
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

module grill(r, R, N) {
    for (k = [0:3]) {
        for (i = [0:N]) {
            rotate([0, 0, k*120 + (i/N - 0.5) * grill_sector_angle])
                translate([R, 0, -1]) cylinder(r=r, h=injection_end_h+2);
        }
    }
}

module grills() {
    grill(grill_hole_r, 47, 16);
    grill(grill_hole_r, 43, 14);
    grill(grill_hole_r, 39, 12);
}

module extraction_cyl() {
    translate([0, 0, cylinder_outer_length - extraction_cyl_h]) difference() {
        cylinder(r=extraction_cyl_r, h=extraction_cyl_h);
        translate([0, 0, -1]) cylinder(r=extraction_hole_r, h=extraction_cyl_offset+2);
    }
}

module extraction_cyl_remove() {
    translate([0, 0, cylinder_outer_length - extraction_cyl_h + extraction_cyl_offset])
    cylinder(r=extraction_cyl_r-extraction_cyl_offset, h=extraction_cyl_h + 1);
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
    translate([0, 0, injection_end_h]) difference() {
        union() {
            translate([0, 0, injection_plug_plate_h - injection_plug_plate_d])
                cylinder(r=injection_plug_plate_r, h=injection_plug_plate_d);
            cylinder(r=injection_plug_bottom_r, h=injection_plug_bottom_h);
            cylinder(r1=injection_plug_r1, r2=injection_plug_r2, h=injection_plug_h);
            translate([0, 0, injection_plug_h])
                cylinder(r=injection_plug_button_r, h=injection_plug_button_h);
        }

        minkowski() {
            waveguides();
            cylinder(r=injection_plug_plate_clearance);
        }
    }
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

module waveguides() {
    for (i = [0:3]) {
        rotate([0, 0, i*120 + 60]) waveguide();
    }
}

module waveguide_holes() {
    for (i = [0:3]) {
        rotate([0, 0, i*120 + 60]) waveguide_hole();
    }
}

module assembly(preview) {
    difference() {
        union() {
            body();
            extraction_cyl();
            waveguides();
            injection_tube();
            injection_plug();
        }

        waveguide_holes();
        injection_tube_hole();
        extraction_cyl_remove();

        if (!preview) {
            grills();
        }
    }
}

module sheath() {
    translate([0, 0, -1])
        cylinder(r = cylinder_r_outer + 1, h = cylinder_outer_length + 2);
}

module circle_surface(r, n) {
    verts = concat([[0, 0, 0]],
        [ for (i = [0:n]) [r*cos(i*360/n), r*sin(i*360/n), 0]]);
    faces = concat([[0, n, 1]], [ for (i = [1:n-1]) [0, i, i+1]]);
    polyhedron(verts, faces);
}

module injection_surface() {
    rotate([0, 0, 60]) translate([injection_tube_R, 0, injection_tube_h])
        circle_surface(injection_tube_r, 20);
}

module extraction_surface() {
    translate([0, 0, cylinder_outer_length - extraction_cyl_h + 1])
        circle_surface(extraction_hole_r + 1, 20);
}

//assembly(true);
//projection(cut=true) translate([0, 0, -100]) assembly(true);

translate([0, 0, -203.8 - injection_end_h - injection_plug_h_total]) {
    if (PARTNO == 1) {
        assembly(false);
    } else if (PARTNO == 2) {
        sheath();
    } else if (PARTNO == 3) {
        injection_surface();
    } else if (PARTNO == 4) {
        extraction_surface();
    } else if (PARTNO == 5) {
        projection(cut=true) translate([0, 0, -100]) assembly(false);
    } else if (PARTNO == 6) {
        projection(cut=true) rotate([0, 0, 0]) assembly(false);
    } else if (PARTNO == 7) {
        projection(cut=true) rotate([90, 0, 0]) assembly(false);
    } else {
         assembly(false);
        //color([1, 0, 0, 0.5]) injection_surface();
        //color([1, 0, 0, 0.5]) extraction_surface();
    }
}

