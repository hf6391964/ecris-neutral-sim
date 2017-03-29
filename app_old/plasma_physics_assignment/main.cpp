#include <iostream>
#include <fstream>
#include <string>

#include "cgal_and_typedefs.h"
#include "electronmodel.h"

int main() {
    const double z1 = -160e-3,
                 z2 = 170e-3,
                 r0 = 39e-3,
                 B0 = 1.07,
                 dt = 5e-12,
                 confinementTime = 1e-6,
                 gridSize = 0.05;
    double speed_10keV = std::sqrt(2.0 * 10e3 / ELECTRON_MASS_EV) *
        SPEED_OF_LIGHT;

    // B0, r0, dt, z1, z2
    ElectronModel model(B0, r0, dt, z1, z2, gridSize, confinementTime);

    // Plot solenoid Bz(r=0, z)
    int nPoints = 100;
    std::ofstream fout("output/solenoid_B_onaxis.csv");
    for (int i = 0; i <= nPoints; i++) {
        double z = z1 + (z2 - z1) * (double)i/nPoints;
        Vector B = ElectronModel::solenoidBfield(0.0, 0.0, z);
        fout << z << CSV_SEP << B.z() << std::endl;
    }
    fout.close();

    // Plot hexapole Bxy(x, y)
    fout.open("output/hexapole_B.csv");
    nPoints = 15;
    for (int i = 0; i <= nPoints; i++) {
        double y = -r0 + 2.0*r0 * (double)i/nPoints;
        for (int j = 0; j <= nPoints; j++) {
            double x = -r0 + 2.0*r0 * (double)j/nPoints;
            Vector B = ElectronModel::hexapoleBfield(x, y, B0, r0);
            fout << x << CSV_SEP << y << CSV_SEP << B.x() << CSV_SEP <<
                B.y() << std::endl;
        }
    }
    fout.close();

    // Total fields Bxz and Byz
    fout.open("output/total_B_xz.csv");
    std::ofstream fout2("output/total_B_yz.csv");
    int nPointsZ = 20;
    int nPointsXY = 8;
    for (int i = 0; i <= nPointsZ; i++) {
        double z = z1 + (z2 - z1) * (double)i/nPointsZ;
        for (int j = 0; j <= nPointsXY; j++) {
            double xy = -r0 + 2.0*r0 * (double)j/nPointsXY;

            Vector Bxz = ElectronModel::totalBfield(Vector(xy, 0.0, z), B0, r0);
            Vector Byz = ElectronModel::totalBfield(Vector(0.0, xy, z), B0, r0);

            fout << xy << CSV_SEP << z << CSV_SEP << Bxz.x() << CSV_SEP <<
                Bxz.z() << std::endl;
            fout2 << xy << CSV_SEP << z << CSV_SEP << Byz.y() << CSV_SEP <<
                Byz.z() << std::endl;
        }
    }
    fout.close();
    fout2.close();

    Rng rng(2016);

    // Track a single particle trajectory to observe stability for a given
    // timestep length
    double simTime = 1e-8;
    const double BnormMax = 0.5;
    double dt_ = 1e-10;
    ElectronModel model1(B0, r0, dt_, z1, z2, gridSize, confinementTime);
    model1.newParticle(speed_10keV,
        [BnormMax](Vector B) { return B.squared_length() < BnormMax*BnormMax; },
        rng);
    fout.open("output/particle_trajectory_10.csv");
    for (double t = 0.0; t < simTime; t += dt_) {
        Vector pos = model1.position();
        fout << t << CSV_SEP << pos.x() << CSV_SEP << pos.y() << CSV_SEP <<
            pos.z() << CSV_SEP << model1.energy() << std::endl;
        model1.moveParticle();
    }
    fout.close();
    dt_ = 1e-11;
    ElectronModel model2(B0, r0, dt_, z1, z2, gridSize, confinementTime);
    rng.seed(2016);
    model2.newParticle(speed_10keV,
        [BnormMax](Vector B) { return B.squared_length() < BnormMax*BnormMax; },
        rng);
    fout.open("output/particle_trajectory_11.csv");
    for (double t = 0.0; t < simTime; t += dt_) {
        Vector pos = model2.position();
        fout << t << CSV_SEP << pos.x() << CSV_SEP << pos.y() << CSV_SEP <<
            pos.z() << CSV_SEP << model2.energy() << std::endl;
        model2.moveParticle();
    }
    fout.close();
    dt_ = 1e-12;
    rng.seed(2016);
    ElectronModel model3(B0, r0, dt_, z1, z2, gridSize, confinementTime);
    model3.newParticle(10e3,
        [BnormMax](Vector B) { return B.squared_length() < BnormMax*BnormMax; },
        rng);
    fout.open("output/particle_trajectory_12.csv");
    for (double t = 0.0; t < simTime; t += dt_) {
        Vector pos = model3.position();
        fout << t << CSV_SEP << pos.x() << CSV_SEP << pos.y() << CSV_SEP <<
            pos.z() << CSV_SEP << model3.energy() << std::endl;
        model3.moveParticle();
    }
    fout.close();

    model.resetCounters();
    // Track multiple particles and capture some statistics of the motion
    unsigned long nParticles = 100000;
    model.runMonoenergeticSimulation(nParticles, 10e3, BnormMax, rng);
    fout.open("output/z1_collision_points.csv");
    for (Vector pos : model.z1CollisionPoints()) {
        fout << pos.x() << CSV_SEP << pos.y() << CSV_SEP << pos.z() << std::endl;
    }
    fout.close();
    fout.open("output/z2_collision_points.csv");
    for (Vector pos : model.z2CollisionPoints()) {
        fout << pos.x() << CSV_SEP << pos.y() << CSV_SEP << pos.z() << std::endl;
    }
    fout.close();
    fout.open("output/cylinder_collision_points.csv");
    for (Vector pos : model.cylinderCollisionPoints()) {
        fout << pos.x() << CSV_SEP << pos.y() << CSV_SEP << pos.z() << std::endl;
    }
    fout.close();

    fout.open("output/nonlost_velocity.csv");
    for (std::tuple<double, double> v : model.nonLostVelocities()) {
        double vr, vt;
        std::tie(vr, vt) = v;
        fout << vr << CSV_SEP << vt << std::endl;
    }
    fout.close();

    fout.open("output/lost_velocity.csv");
    for (std::tuple<double, double> v : model.lostVelocities()) {
        double vr, vt;
        std::tie(vr, vt) = v;
        fout << vr << CSV_SEP << vt << std::endl;
    }
    fout.close();

    fout.open("output/energy.txt");
    fout << "Mean energy: " << model.meanEnergy() <<
        " eV, standard deviation: " << model.energyStdDev() <<
        " eV" << std::endl;
    fout.close();

    return 0;
}

