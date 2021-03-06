/********************************************************************************************************************************************
 * Saras
 * 
 * Copyright (C) 2019, Mahendra K. Verma
 *
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     1. Redistributions of source code must retain the above copyright
 *        notice, this list of conditions and the following disclaimer.
 *     2. Redistributions in binary form must reproduce the above copyright
 *        notice, this list of conditions and the following disclaimer in the
 *        documentation and/or other materials provided with the distribution.
 *     3. Neither the name of the copyright holder nor the
 *        names of its contributors may be used to endorse or promote products
 *        derived from this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 ********************************************************************************************************************************************
 */
/*! \file scalar_d3.cc
 *
 *  \brief Definitions of functions for 3D computations with the scalar class.
 *  \sa scalar.h
 *  \author Shashwat Bhattacharya, Roshan Samuel
 *  \date Nov 2019
 *  \copyright New BSD License
 *
 ********************************************************************************************************************************************
 */

#include <sys/time.h>
#include <ctime>

#include "scalar.h"

/**
 ********************************************************************************************************************************************
 * \brief   Constructor of the scalar_d3 class derived from the base scalar class
 *
 *          The constructor passes its arguments to the base scalar class and then initializes all the scalar and
 *          vector fields necessary for solving the NS equations.
 *          The various coefficients for solving the equations are also set by a call to the \ref setCoefficients function.
 *          Based on the problem type specified by the user in the parameters file, and stored by the \ref parser class as
 *          \ref parser#probType "probType", the appropriate boundary conditions are specified.
 *
 * \param   mesh is a const reference to the global data contained in the grid class
 * \param   solParam is a const reference to the user-set parameters contained in the parser class
 ********************************************************************************************************************************************
 */
scalar_d3::scalar_d3(const grid &mesh, const parser &solParam, parallel &mpiParam):
            scalar(mesh, solParam, mpiParam),
            mgSolver(mesh, inputParams)
{
    // SET VALUES OF COEFFICIENTS USED FOR COMPUTING LAPLACIAN
    setCoefficients();

    // INITIALIZE VARIABLES
    if (inputParams.restartFlag) {
        // Fields to be read from HDF5 file are passed to reader class as a vector
        std::vector<field> readFields;

        // Populate the vector with required fields
        readFields.push_back(V.Vx);
        readFields.push_back(V.Vy);
        readFields.push_back(V.Vz);
        readFields.push_back(P.F);
        readFields.push_back(T.F);

        // Initialize reader object
        reader dataReader(mesh, readFields);

        time = dataReader.readData();

        // Abort if this time is greater than the final time specified by the user
        if (time >= inputParams.tMax) {
            if (mesh.rankData.rank == 0) {
                std::cout << "ERROR: Restart file is starting from a point beyond the final time specified. Aborting" << std::endl;
            }
            MPI_Finalize();
            exit(0);
        }
    } else {
        P = 1.0;
        V = 0.0;
     for (int i=T.F.fBulk.lbound(0); i <= T.F.fBulk.ubound(0); i++) {
        for (int j=T.F.fBulk.lbound(1); j <= T.F.fBulk.ubound(1); j++) {
            for (int k=T.F.fBulk.lbound(2); k <= T.F.fBulk.ubound(2); k++) {
                T.F.F(i, j, k) = (1.0 - mesh.zStaggr(k)/mesh.zLen);

            }
        }
    }
        
        omega_x = 0.0, omega_y = 0.0, omega_z = 0.0;
        instomega_x = 0.0, instomega_y = 0.0, instomega_z = 0.0;

        Hlctx = 0.0, Hlcty = 0.0, Hlctz = 0.0;      
        instHlctx = 0.0, instHlcty = 0.0, instHlctz = 0.0;      

        time = 0.0;
    }

    checkPeriodic();

    // Initialize velocity and temperature boundary conditions
    initVBCs();
    initTBCs();

    // Initialize velocity and temperature forcing fields
    initVForcing();
    initTForcing();

    // Impose boundary conditions on velocity and temperature fields
    V.imposeBCs();
    T.imposeBCs();
}


void scalar_d3::solvePDE() {
#ifdef TIME_RUN
    visc_time = 0.0;
    nlin_time = 0.0;
    intr_time = 0.0;
    impl_time = 0.0;
    prhs_time = 0.0;
    pois_time = 0.0;

#else
    real fwTime, prTime, rsTime;

    // Set dt equal to input time step
    dt = inputParams.tStp;

    // Fields to be written into HDF5 file are passed to writer class as a vector
    std::vector<field> writeFields;

    // Populate the vector with required fields
    writeFields.push_back(V.Vx);
    writeFields.push_back(V.Vy);
    writeFields.push_back(V.Vz);
    writeFields.push_back(P.F);
    writeFields.push_back(T.F);
    //writeFields.push_back(instHlctx);
    //writeFields.push_back(instHlcty);
    //writeFields.push_back(instHlctz);
    //writeFields.push_back(Helicity.F);

    // Initialize writer object
    writer dataWriter(mesh, writeFields);

    // Initialize probes
    if (inputParams.readProbes) {
        dataProbe = new probes(mesh, writeFields);
    }

    // Output file containing the time series of various variables
    tseries tsWriter(mesh, V, P, time, dt);

    // FILE WRITING TIME
    fwTime = time;

    // FIELD PROBING TIME
    prTime = time;

    // RESTART FILE WRITING TIME
    rsTime = time;
    

    //#################################time-avg-parameters########################

	//if (inputParams.AvgSolnFlag) {
//params 
    int avgCounter = 0;
    //double acStrt = 1.0;


// Averaged Fields to be written into HDF5 file are passed to writer class as a vector
    std::vector<field> writeAvgFields;

    // Populate the vector with required scalar fields
    writeAvgFields.push_back(omega_x.F);
    writeAvgFields.push_back(omega_y.F);
    writeAvgFields.push_back(omega_z.F);
    writeAvgFields.push_back(Hlctx.F);
    writeAvgFields.push_back(Hlcty.F);
    writeAvgFields.push_back(Hlctz.F);

    // Initialize writer object
    writer avgDataWriter(mesh, writeAvgFields);    
    
	//}
    
#endif

    timeStepCount = 0;

#ifndef TIME_RUN
    // COMPUTE ENERGY AND DIVERGENCE FOR THE INITIAL CONDITION
    tsWriter.writeTSData(T, nu, kappa);

    // WRITE DATA AT t = 0 OR INCREMENT INTERVAL IF RESTARTING
    if (inputParams.restartFlag) {
        int tCount, fCount;

        tCount = int(time/inputParams.tStp);

        fCount = int(inputParams.fwInt/inputParams.tStp);
        fwTime = roundNum(tCount, fCount)*inputParams.tStp;

        fCount = int(inputParams.prInt/inputParams.tStp);
        prTime = roundNum(tCount, fCount)*inputParams.tStp;

        fCount = int(inputParams.rsInt/inputParams.tStp);
        rsTime = roundNum(tCount, fCount)*inputParams.tStp;
    }

    switch (inputParams.solnFormat) {
        case 1: dataWriter.writeSolution(time);
            break;
        case 2: dataWriter.writeTarang(time);
            break;
        default: dataWriter.writeSolution(time);
    }
    fwTime += inputParams.fwInt;

    if (inputParams.readProbes) {
        dataProbe->probeData(time);
        prTime += inputParams.prInt;
    }

    rsTime += inputParams.rsInt;
#endif

    // TIME-INTEGRATION LOOP
    while (true) {
        // MAIN FUNCTION CALLED IN EACH LOOP TO UPDATE THE FIELDS AT EACH TIME-STEP
        timeAdvance();

        if (inputParams.useCFL) {
            V.computeTStp(dt);
            if (dt > inputParams.tStp) {
                dt = inputParams.tStp;
            }
        }

        timeStepCount += 1;
        time += dt;

#ifndef TIME_RUN
        if (timeStepCount % inputParams.ioCnt == 0) {
            tsWriter.writeTSData(T, nu, kappa);
        }

        if (inputParams.readProbes and std::abs(prTime - time) < 0.5*dt) {
            dataProbe->probeData(time);
            prTime += inputParams.prInt;
        }

        if (std::abs(fwTime - time) < 0.5*dt) {
            switch (inputParams.solnFormat) {
                case 1: dataWriter.writeSolution(time);
                    break;
                case 2: dataWriter.writeTarang(time);
                    break;
                default: dataWriter.writeSolution(time);
            }
            fwTime += inputParams.fwInt;
        }

	if (inputParams.AvgSolnFlag) {
       if (time >= inputParams.AvgStartTime) {
            avgCounter += 1;
            omega_x += instomega_x;
            omega_y += instomega_y;
            omega_z += instomega_z;
            Hlctx += instHlctx;
            Hlcty += instHlcty;
            Hlctz += instHlctz;

            
			//std::cout << awTime << std::endl;            
            //if (std::abs(awTime - time) < 0.5*dt) {
	    if (std::abs(inputParams.AvgWriteTime - time) < 0.5*dt) {					
                 omega_x *= 1.0/avgCounter;
                 omega_y *= 1.0/avgCounter; 
                 omega_z *= 1.0/avgCounter; 
                 Hlctx *= 1.0/avgCounter;
                 Hlcty *= 1.0/avgCounter; 
                 Hlctz *= 1.0/avgCounter;                  
                 
                 avgDataWriter.writeAvgSolution(time);
            }                    
        }    
	}


        if (std::abs(rsTime - time) < 0.5*dt) {
            dataWriter.writeRestart(time);
            rsTime += inputParams.rsInt;
        }
#endif

        if (std::abs(inputParams.tMax - time) < 0.5*dt) {
            break;
        }
    }

    // WRITE THE OUTPUT OF THE TIMING RUN
#ifdef TIME_RUN
    if (mesh.rankData.rank == 0) {
        std::cout << std::left << std::setw(50) << "Time taken in computing viscous terms: "            << std::fixed << std::setprecision(6) << visc_time << std::endl;
        std::cout << std::left << std::setw(50) << "Time taken in computing non-linear terms: "         << std::fixed << std::setprecision(6) << nlin_time << std::endl;
        std::cout << std::left << std::setw(50) << "Time taken in computing intermediate steps: "       << std::fixed << std::setprecision(6) << intr_time << std::endl;
        std::cout << std::left << std::setw(50) << "Time taken in computing velocities implicitly: "    << std::fixed << std::setprecision(6) << impl_time << std::endl;
        std::cout << std::left << std::setw(50) << "Time taken in computing RHS for poisson solver: "   << std::fixed << std::setprecision(6) << prhs_time << std::endl;
        std::cout << std::left << std::setw(50) << "Time taken by poisson solver: "                     << std::fixed << std::setprecision(6) << pois_time << std::endl;
    }
#endif
}


void scalar_d3::timeAdvance() {
#ifdef TIME_RUN
    struct timeval begin, end;
#endif

    nseRHS = 0.0;
    tmpRHS = 0.0;

    // First compute the explicit part of the semi-implicit viscous term and divide it by Re
#ifdef TIME_RUN
    gettimeofday(&begin, NULL);
    V.computeDiff(nseRHS);
    nseRHS *= nu;
    gettimeofday(&end, NULL);
    visc_time += ((end.tv_sec - begin.tv_sec)*1000000u + end.tv_usec - begin.tv_usec)/1.e6;
#else
    V.computeDiff(nseRHS);
    nseRHS *= nu;
#endif

    // Compute the non-linear term and subtract it from the RHS
#ifndef TIME_RUN
    V.computeNLin(V, nseRHS);
#else
    gettimeofday(&begin, NULL);
    V.computeNLin(V, nseRHS);
    gettimeofday(&end, NULL);
    nlin_time += ((end.tv_sec - begin.tv_sec)*1000000u + end.tv_usec - begin.tv_usec)/1.e6;

    gettimeofday(&begin, NULL);
#endif

    // Add the velocity forcing term
    V.vForcing->addForcing(nseRHS);

    // Subtract the pressure gradient term
    pressureGradient = 0.0;
    P.gradient(pressureGradient, V);
    nseRHS -= pressureGradient;

    // Multiply the entire RHS with dt and add the velocity of previous time-step to advance by explicit Euler method
    nseRHS *= dt;
    nseRHS += V;

#ifdef TIME_RUN
    gettimeofday(&end, NULL);
    intr_time += ((end.tv_sec - begin.tv_sec)*1000000u + end.tv_usec - begin.tv_usec)/1.e6;
#endif

    // Synchronize the RHS term across all processors by updating its sub-domain pads
    nseRHS.syncData();

    // Using the RHS term computed, compute the guessed velocity of CN method iteratively (and store it in V)
#ifdef TIME_RUN
    gettimeofday(&begin, NULL);
#endif
    solveVx();
    solveVy();
    solveVz();

#ifdef TIME_RUN
    gettimeofday(&end, NULL);
    impl_time += ((end.tv_sec - begin.tv_sec)*1000000u + end.tv_usec - begin.tv_usec)/1.e6;
#endif

    // Calculate the rhs for the poisson solver (mgRHS) using the divergence of guessed velocity in V
#ifdef TIME_RUN
    gettimeofday(&begin, NULL);
    V.divergence(mgRHS, P);
    mgRHS *= 1.0/dt;
    gettimeofday(&end, NULL);
    prhs_time += ((end.tv_sec - begin.tv_sec)*1000000u + end.tv_usec - begin.tv_usec)/1.e6;
#else
    V.divergence(mgRHS, P);
    mgRHS *= 1.0/dt;
#endif

    // Using the calculated mgRHS, evaluate pressure correction (Pp) using multi-grid method
#ifdef TIME_RUN
    gettimeofday(&begin, NULL);
    mgSolver.mgSolve(Pp, mgRHS);
    gettimeofday(&end, NULL);
    pois_time += ((end.tv_sec - begin.tv_sec)*1000000u + end.tv_usec - begin.tv_usec)/1.e6;
#else
    mgSolver.mgSolve(Pp, mgRHS);
#endif

    // Synchronise the pressure correction term across processors
    Pp.syncData();

    // Add the pressure correction term to the pressure field of previous time-step, P
    P += Pp;

    // Finally get the velocity field at end of time-step by subtracting the gradient of pressure correction from V
    Pp.gradient(pressureGradient, V);
    pressureGradient *= dt;
    V -= pressureGradient;

    // Next compute the explicit part of the semi-implicit diffusion term and multiply it by kappa
    T.computeDiff(tmpRHS);
    tmpRHS *= kappa;

    // COMPUTE THE CONVECTIVE DERIVATIVE AND SUBTRACT IT FROM THE CALCULATED DIFFUSION TERMS OF RHS IN tmpRHS
    T.computeNLin(V, tmpRHS);

    // Add the scalar forcing term
    T.tForcing->addForcing(tmpRHS);

    // Multiply the entire RHS with dt and add the temperature of previous time-step to advance by explicit Euler method
    tmpRHS *= dt;
    tmpRHS += T;

    // Synchronize the RHS term across all processors by updating its sub-domain pads
    tmpRHS.syncData();

    // Using the RHS term computed, compute the guessed temperature of CN method iteratively (and store it in T)
    solveT();

    FindHelicity();

    // Impose boundary conditions on the updated velocity field, V
    V.imposeBCs();

    // Impose boundary conditions on the updated temperature field, T
    T.imposeBCs();
}




void scalar_d3::solveVx() {
    int iterCount = 0;
    real maxError = 0.0;

    while (true) {
#pragma omp parallel for num_threads(inputParams.nThreads) default(none)
        for (int iX = V.Vx.fBulk.lbound(0); iX <= V.Vx.fBulk.ubound(0); iX++) {
            for (int iY = V.Vx.fBulk.lbound(1); iY <= V.Vx.fBulk.ubound(1); iY++) {
                for (int iZ = V.Vx.fBulk.lbound(2); iZ <= V.Vx.fBulk.ubound(2); iZ++) {
                    guessedVelocity.Vx(iX, iY, iZ) = ((hy2hz2 * mesh.xix2Colloc(iX) * (V.Vx.F(iX+1, iY, iZ) + V.Vx.F(iX-1, iY, iZ)) +
                                                       hz2hx2 * mesh.ety2Staggr(iY) * (V.Vx.F(iX, iY+1, iZ) + V.Vx.F(iX, iY-1, iZ)) +
                                                       hx2hy2 * mesh.ztz2Staggr(iZ) * (V.Vx.F(iX, iY, iZ+1) + V.Vx.F(iX, iY, iZ-1))) *
                            dt * nu / ( hx2hy2hz2 * 2.0) + nseRHS.Vx(iX, iY, iZ)) /
                     (1.0 + dt * nu * ((hy2hz2 * mesh.xix2Colloc(iX) + hz2hx2 * mesh.ety2Staggr(iY) + hx2hy2 * mesh.ztz2Staggr(iZ)))/hx2hy2hz2);
                }
            }
        }

        V.Vx.F = guessedVelocity.Vx;

        V.imposeVxBC();

#pragma omp parallel for num_threads(inputParams.nThreads) default(none)
        for (int iX = V.Vx.fBulk.lbound(0); iX <= V.Vx.fBulk.ubound(0); iX++) {
            for (int iY = V.Vx.fBulk.lbound(1); iY <= V.Vx.fBulk.ubound(1); iY++) {
                for (int iZ = V.Vx.fBulk.lbound(2); iZ <= V.Vx.fBulk.ubound(2); iZ++) {
                    velocityLaplacian.Vx(iX, iY, iZ) = V.Vx.F(iX, iY, iZ) - 0.5 * dt * nu * (
                              mesh.xix2Colloc(iX) * (V.Vx.F(iX+1, iY, iZ) - 2.0 * V.Vx.F(iX, iY, iZ) + V.Vx.F(iX-1, iY, iZ)) / (hx * hx) +
                              mesh.ety2Staggr(iY) * (V.Vx.F(iX, iY+1, iZ) - 2.0 * V.Vx.F(iX, iY, iZ) + V.Vx.F(iX, iY-1, iZ)) / (hy * hy) +
                              mesh.ztz2Staggr(iZ) * (V.Vx.F(iX, iY, iZ+1) - 2.0 * V.Vx.F(iX, iY, iZ) + V.Vx.F(iX, iY, iZ-1)) / (hz * hz));
                }
            }
        }

        velocityLaplacian.Vx(V.Vx.fBulk) = abs(velocityLaplacian.Vx(V.Vx.fBulk) - nseRHS.Vx(V.Vx.fBulk));

        maxError = velocityLaplacian.vxMax();

        if (maxError < inputParams.cnTolerance) {
            break;
        }

        iterCount += 1;

        if (iterCount > maxIterations) {
            if (mesh.rankData.rank == 0) {
                std::cout << "ERROR: Jacobi iterations for solution of Vx not converging. Aborting" << std::endl;
            }
            MPI_Finalize();
            exit(0);
        }
    }
}


void scalar_d3::solveVy() {
    int iterCount = 0;
    real maxError = 0.0;

    while (true) {
#pragma omp parallel for num_threads(inputParams.nThreads) default(none)
        for (int iX = V.Vy.fBulk.lbound(0); iX <= V.Vy.fBulk.ubound(0); iX++) {
            for (int iY = V.Vy.fBulk.lbound(1); iY <= V.Vy.fBulk.ubound(1); iY++) {
                for (int iZ = V.Vy.fBulk.lbound(2); iZ <= V.Vy.fBulk.ubound(2); iZ++) {
                    guessedVelocity.Vy(iX, iY, iZ) = ((hy2hz2 * mesh.xix2Staggr(iX) * (V.Vy.F(iX+1, iY, iZ) + V.Vy.F(iX-1, iY, iZ)) +
                                                       hz2hx2 * mesh.ety2Colloc(iY) * (V.Vy.F(iX, iY+1, iZ) + V.Vy.F(iX, iY-1, iZ)) +
                                                       hx2hy2 * mesh.ztz2Staggr(iZ) * (V.Vy.F(iX, iY, iZ+1) + V.Vy.F(iX, iY, iZ-1))) *
                            dt * nu / ( hx2hy2hz2 * 2.0) + nseRHS.Vy(iX, iY, iZ)) /
                     (1.0 + dt * nu * ((hy2hz2 * mesh.xix2Staggr(iX) + hz2hx2 * mesh.ety2Colloc(iY) + hx2hy2 * mesh.ztz2Staggr(iZ)))/hx2hy2hz2);
                }
            }
        }

        V.Vy.F = guessedVelocity.Vy;

        V.imposeVyBC();

#pragma omp parallel for num_threads(inputParams.nThreads) default(none)
        for (int iX = V.Vy.fBulk.lbound(0); iX <= V.Vy.fBulk.ubound(0); iX++) {
            for (int iY = V.Vy.fBulk.lbound(1); iY <= V.Vy.fBulk.ubound(1); iY++) {
                for (int iZ = V.Vy.fBulk.lbound(2); iZ <= V.Vy.fBulk.ubound(2); iZ++) {
                    velocityLaplacian.Vy(iX, iY, iZ) = V.Vy.F(iX, iY, iZ) - 0.5 * dt * nu * (
                              mesh.xix2Staggr(iX) * (V.Vy.F(iX+1, iY, iZ) - 2.0 * V.Vy.F(iX, iY, iZ) + V.Vy.F(iX-1, iY, iZ)) / (hx * hx) +
                              mesh.ety2Colloc(iY) * (V.Vy.F(iX, iY+1, iZ) - 2.0 * V.Vy.F(iX, iY, iZ) + V.Vy.F(iX, iY-1, iZ)) / (hy * hy) +
                              mesh.ztz2Staggr(iZ) * (V.Vy.F(iX, iY, iZ+1) - 2.0 * V.Vy.F(iX, iY, iZ) + V.Vy.F(iX, iY, iZ-1)) / (hz * hz));
                }
            }
        }

        velocityLaplacian.Vy(V.Vy.fBulk) = abs(velocityLaplacian.Vy(V.Vy.fBulk) - nseRHS.Vy(V.Vy.fBulk));

        maxError = velocityLaplacian.vyMax();

        if (maxError < inputParams.cnTolerance) {
            break;
        }

        iterCount += 1;

        if (iterCount > maxIterations) {
            if (mesh.rankData.rank == 0) {
                std::cout << "ERROR: Jacobi iterations for solution of Vy not converging. Aborting" << std::endl;
            }
            MPI_Finalize();
            exit(0);
        }
    }
}


void scalar_d3::solveVz() {
    int iterCount = 0;
    real maxError = 0.0;

    while (true) {
#pragma omp parallel for num_threads(inputParams.nThreads) default(none)
        for (int iX = V.Vz.fBulk.lbound(0); iX <= V.Vz.fBulk.ubound(0); iX++) {
            for (int iY = V.Vz.fBulk.lbound(1); iY <= V.Vz.fBulk.ubound(1); iY++) {
                for (int iZ = V.Vz.fBulk.lbound(2); iZ <= V.Vz.fBulk.ubound(2); iZ++) {
                    guessedVelocity.Vz(iX, iY, iZ) = ((hy2hz2 * mesh.xix2Staggr(iX) * (V.Vz.F(iX+1, iY, iZ) + V.Vz.F(iX-1, iY, iZ)) +
                                                       hz2hx2 * mesh.ety2Staggr(iY) * (V.Vz.F(iX, iY+1, iZ) + V.Vz.F(iX, iY-1, iZ)) +
                                                       hx2hy2 * mesh.ztz2Colloc(iZ) * (V.Vz.F(iX, iY, iZ+1) + V.Vz.F(iX, iY, iZ-1))) *
                            dt * nu / ( hx2hy2hz2 * 2.0) + nseRHS.Vz(iX, iY, iZ)) /
                     (1.0 + dt * nu * ((hy2hz2 * mesh.xix2Staggr(iX) + hz2hx2 * mesh.ety2Staggr(iY) + hx2hy2 * mesh.ztz2Colloc(iZ)))/hx2hy2hz2);
                }
            }
        }

        V.Vz.F = guessedVelocity.Vz;

        V.imposeVzBC();

#pragma omp parallel for num_threads(inputParams.nThreads) default(none)
        for (int iX = V.Vz.fBulk.lbound(0); iX <= V.Vz.fBulk.ubound(0); iX++) {
            for (int iY = V.Vz.fBulk.lbound(1); iY <= V.Vz.fBulk.ubound(1); iY++) {
                for (int iZ = V.Vz.fBulk.lbound(2); iZ <= V.Vz.fBulk.ubound(2); iZ++) {
                    velocityLaplacian.Vz(iX, iY, iZ) = V.Vz.F(iX, iY, iZ) - 0.5 * dt * nu * (
                              mesh.xix2Staggr(iX) * (V.Vz.F(iX+1, iY, iZ) - 2.0 * V.Vz.F(iX, iY, iZ) + V.Vz.F(iX-1, iY, iZ)) / (hx * hx) +
                              mesh.ety2Staggr(iY) * (V.Vz.F(iX, iY+1, iZ) - 2.0 * V.Vz.F(iX, iY, iZ) + V.Vz.F(iX, iY-1, iZ)) / (hy * hy) +
                              mesh.ztz2Colloc(iZ) * (V.Vz.F(iX, iY, iZ+1) - 2.0 * V.Vz.F(iX, iY, iZ) + V.Vz.F(iX, iY, iZ-1)) / (hz * hz));
                }
            }
        }

        velocityLaplacian.Vz(V.Vz.fBulk) = abs(velocityLaplacian.Vz(V.Vz.fBulk) - nseRHS.Vz(V.Vz.fBulk));

        maxError = velocityLaplacian.vzMax();

        if (maxError < inputParams.cnTolerance) {
            break;
        }

        iterCount += 1;

        if (iterCount > maxIterations) {
            if (mesh.rankData.rank == 0) {
                std::cout << "ERROR: Jacobi iterations for solution of Vz not converging. Aborting" << std::endl;
            }
            MPI_Finalize();
            exit(0);
        }
    }
}


void scalar_d3::solveT() {
    int iterCount = 0;
    real maxError = 0.0;

    while (true) {
#pragma omp parallel for num_threads(inputParams.nThreads) default(none)
        for (int iX = T.F.fBulk.lbound(0); iX <= T.F.fBulk.ubound(0); iX++) {
            for (int iY = T.F.fBulk.lbound(1); iY <= T.F.fBulk.ubound(1); iY++) {
                for (int iZ = T.F.fBulk.lbound(2); iZ <= T.F.fBulk.ubound(2); iZ++) {
                    guessedScalar.F(iX, iY, iZ) = ((hy2hz2 * mesh.xix2Staggr(iX) * (T.F.F(iX+1, iY, iZ) + T.F.F(iX-1, iY, iZ)) +
                                                    hz2hx2 * mesh.ety2Staggr(iY) * (T.F.F(iX, iY+1, iZ) + T.F.F(iX, iY-1, iZ)) +
                                                    hx2hy2 * mesh.ztz2Colloc(iZ) * (T.F.F(iX, iY, iZ+1) + T.F.F(iX, iY, iZ-1))) *
                             dt * kappa / ( hx2hy2hz2 * 2.0) + tmpRHS.F(iX, iY, iZ)) /
                      (1.0 + dt * kappa * ((hy2hz2 * mesh.xix2Staggr(iX) + hz2hx2 * mesh.ety2Staggr(iY) + hx2hy2 * mesh.ztz2Colloc(iZ)))/hx2hy2hz2);
                }
            }
        }

        T = guessedScalar;

        T.imposeBCs();

#pragma omp parallel for num_threads(inputParams.nThreads) default(none)
        for (int iX = T.F.fBulk.lbound(0); iX <= T.F.fBulk.ubound(0); iX++) {
            for (int iY = T.F.fBulk.lbound(1); iY <= T.F.fBulk.ubound(1); iY++) {
                for (int iZ = T.F.fBulk.lbound(2); iZ <= T.F.fBulk.ubound(2); iZ++) {
                    scalarLaplacian.F(iX, iY, iZ) = T.F.F(iX, iY, iZ) - 0.5 * dt * kappa * (
                           mesh.xix2Staggr(iX) * (T.F.F(iX+1, iY, iZ) - 2.0 * T.F.F(iX, iY, iZ) + T.F.F(iX-1, iY, iZ)) / (hx * hx) +
                           mesh.ety2Staggr(iY) * (T.F.F(iX, iY+1, iZ) - 2.0 * T.F.F(iX, iY, iZ) + T.F.F(iX, iY-1, iZ)) / (hy * hy) +
                           mesh.ztz2Colloc(iZ) * (T.F.F(iX, iY, iZ+1) - 2.0 * T.F.F(iX, iY, iZ) + T.F.F(iX, iY, iZ-1)) / (hz * hz));
                }
            }
        }

        scalarLaplacian.F(T.F.fBulk) = abs(scalarLaplacian.F(T.F.fBulk) - tmpRHS.F(T.F.fBulk));

        maxError = scalarLaplacian.fxMax();

        if (maxError < inputParams.cnTolerance) {
            break;
        }

        iterCount += 1;

        if (iterCount > maxIterations) {
            if (mesh.rankData.rank == 0) {
                std::cout << "ERROR: Jacobi iterations for solution of T not converging. Aborting" << std::endl;
            }
            MPI_Finalize();
            exit(0);
        }
    }
}



void scalar_d3::FindHelicity() {
    for (int iX = T.F.fBulk.lbound(0); iX <= T.F.fBulk.ubound(0); iX++) {
        for (int iY = T.F.fBulk.lbound(1); iY <= T.F.fBulk.ubound(1); iY++) {
            for (int iZ = T.F.fBulk.lbound(2); iZ <= T.F.fBulk.ubound(2); iZ++) {
                instomega_x.F(iX, iY, iZ) = ((V.Vz.F(iX, iY+1, iZ) - V.Vz.F(iX, iY-1, iZ))/(2.0*hy)) - ((V.Vy.F(iX, iY, iZ+1) - V.Vy.F(iX, iY, iZ-1))/(2.0*hz));
                instomega_y.F(iX, iY, iZ) = ((V.Vx.F(iX, iY, iZ+1) - V.Vx.F(iX, iY, iZ-1))/(2.0*hz)) - ((V.Vz.F(iX+1, iY, iZ) - V.Vz.F(iX-1, iY, iZ))/(2.0*hx));
                instomega_z.F(iX, iY, iZ) = ((V.Vy.F(iX+1, iY, iZ) - V.Vy.F(iX-1, iY, iZ))/(2.0*hx)) - ((V.Vx.F(iX, iY+1, iZ) - V.Vx.F(iX, iY-1, iZ))/(2.0*hy));
                instHlctx.F(iX, iY, iZ) = V.Vx.F(iX, iY, iZ)*instomega_x.F(iX, iY, iZ); 
                instHlcty.F(iX, iY, iZ) = V.Vy.F(iX, iY, iZ)*instomega_y.F(iX, iY, iZ);
                instHlctz.F(iX, iY, iZ) = V.Vz.F(iX, iY, iZ)*instomega_z.F(iX, iY, iZ);


            }
        }
    }
}




real scalar_d3::testPeriodic() {
    real xCoord = 0.0;
    real yCoord = 0.0;
    real zCoord = 0.0;

    nseRHS = 0.0;
    V = 0.0;

    for (int i=V.Vx.F.lbound(0); i <= V.Vx.F.ubound(0); i++) {
        for (int j=V.Vx.F.lbound(1); j <= V.Vx.F.ubound(1); j++) {
            for (int k=V.Vx.F.lbound(2); k <= V.Vx.F.ubound(2); k++) {
                V.Vx.F(i, j, k) = sin(2.0*M_PI*mesh.xColloc(i)/mesh.xLen)*
                                  cos(2.0*M_PI*mesh.yStaggr(j)/mesh.yLen)*
                                  cos(2.0*M_PI*mesh.zStaggr(k)/mesh.zLen);
                nseRHS.Vx(i, j, k) = V.Vx.F(i, j, k);
            }
        }
    }

    for (int i=V.Vy.F.lbound(0); i <= V.Vy.F.ubound(0); i++) {
        for (int j=V.Vy.F.lbound(1); j <= V.Vy.F.ubound(1); j++) {
            for (int k=V.Vy.F.lbound(2); k <= V.Vy.F.ubound(2); k++) {
                V.Vy.F(i, j, k) = -cos(2.0*M_PI*mesh.xStaggr(i)/mesh.xLen)*
                                   sin(2.0*M_PI*mesh.yColloc(j)/mesh.yLen)*
                                   cos(2.0*M_PI*mesh.zStaggr(k)/mesh.zLen);
                nseRHS.Vy(i, j, k) = V.Vy.F(i, j, k);
            }
        }
    }

    V.Vz.F = 0.0;
    nseRHS.Vz = 0.0;

    // EXPECTED VALUES IN THE PAD REGIONS IF DATA TRANSFER HAPPENS WITH NO HITCH
    // X-VELOCITY IN LEFT AND RIGHT PADS
    for (int iX = 1; iX <= mesh.padWidths(0); iX++) {
        for (int iY = V.Vx.fCore.lbound(1); iY <= V.Vx.fCore.ubound(1); iY += iX) {
            for (int iZ = V.Vx.fCore.lbound(2); iZ <= V.Vx.fCore.ubound(2); iZ += iX) {
                xCoord = mesh.xColloc(V.Vx.fCore.lbound(0)) - (mesh.xColloc(V.Vx.fCore.lbound(0) + iX) - mesh.xColloc(V.Vx.fCore.lbound(0)));
                nseRHS.Vx(V.Vx.fCore.lbound(0) - iX, iY, iZ) = sin(2.0*M_PI*xCoord/mesh.xLen)*
                                                               cos(2.0*M_PI*mesh.yStaggr(iY)/mesh.yLen)*
                                                               cos(2.0*M_PI*mesh.zStaggr(iZ)/mesh.zLen);

                xCoord = mesh.xColloc(V.Vx.fCore.ubound(0)) + (mesh.xColloc(V.Vx.fCore.ubound(0)) - mesh.xColloc(V.Vx.fCore.ubound(0) - iX));
                nseRHS.Vx(V.Vx.fCore.ubound(0) + iX, iY, iZ) = sin(2.0*M_PI*xCoord/mesh.xLen)*
                                                               cos(2.0*M_PI*mesh.yStaggr(iY)/mesh.yLen)*
                                                               cos(2.0*M_PI*mesh.zStaggr(iZ)/mesh.zLen);
            }
        }
    }

    // X-VELOCITY IN FRONT AND BACK PADS
    for (int iY = 1; iY <= mesh.padWidths(1); iY++) {
        for (int iX = V.Vx.fCore.lbound(0); iX <= V.Vx.fCore.ubound(0); iX += iY) {
            for (int iZ = V.Vx.fCore.lbound(2); iZ <= V.Vx.fCore.ubound(2); iZ += iY) {
                yCoord = mesh.yStaggr(V.Vx.fCore.lbound(1)) - (mesh.yStaggr(V.Vx.fCore.lbound(1) + iY) - mesh.yStaggr(V.Vx.fCore.lbound(1)));
                nseRHS.Vx(iX, V.Vx.fCore.lbound(1) - iY, iZ) = sin(2.0*M_PI*mesh.xColloc(iX)/mesh.xLen)*
                                                               cos(2.0*M_PI*yCoord/mesh.yLen)*
                                                               cos(2.0*M_PI*mesh.zStaggr(iZ)/mesh.zLen);

                yCoord = mesh.yStaggr(V.Vx.fCore.ubound(1)) + (mesh.yStaggr(V.Vx.fCore.ubound(1)) - mesh.yStaggr(V.Vx.fCore.ubound(1) - iY));
                nseRHS.Vx(iX, V.Vx.fCore.ubound(1) + iY, iZ) = sin(2.0*M_PI*mesh.xColloc(iX)/mesh.xLen)*
                                                               cos(2.0*M_PI*yCoord/mesh.yLen)*
                                                               cos(2.0*M_PI*mesh.zStaggr(iZ)/mesh.zLen);
            }
        }
    }

    // X-VELOCITY IN BOTTOM AND TOP PADS
    for (int iZ = 1; iZ <= mesh.padWidths(2); iZ++) {
        for (int iX = V.Vx.fCore.lbound(0); iX <= V.Vx.fCore.ubound(0); iX += iZ) {
            for (int iY = V.Vx.fCore.lbound(1); iY <= V.Vx.fCore.ubound(1); iY += iZ) {
                zCoord = mesh.zStaggr(V.Vx.fCore.lbound(2)) - (mesh.zStaggr(V.Vx.fCore.lbound(2) + iZ) - mesh.zStaggr(V.Vx.fCore.lbound(2)));
                nseRHS.Vx(iX, iY, V.Vx.fCore.lbound(2) - iZ) = sin(2.0*M_PI*mesh.xColloc(iX)/mesh.xLen)*
                                                               cos(2.0*M_PI*mesh.yStaggr(iY)/mesh.yLen)*
                                                               cos(2.0*M_PI*zCoord/mesh.zLen);

                zCoord = mesh.zStaggr(V.Vx.fCore.ubound(2)) + (mesh.zStaggr(V.Vx.fCore.ubound(2)) - mesh.zStaggr(V.Vx.fCore.ubound(2) - iZ));
                nseRHS.Vx(iX, iY, V.Vx.fCore.ubound(2) + iZ) = sin(2.0*M_PI*mesh.xColloc(iX)/mesh.xLen)*
                                                               cos(2.0*M_PI*mesh.yStaggr(iY)/mesh.yLen)*
                                                               cos(2.0*M_PI*zCoord/mesh.zLen);
            }
        }
    }

    // Y-VELOCITY IN LEFT AND RIGHT PADS
    for (int iX = 1; iX <= mesh.padWidths(0); iX++) {
        for (int iY = V.Vy.fCore.lbound(1); iY <= V.Vy.fCore.ubound(1); iY += iX) {
            for (int iZ = V.Vy.fCore.lbound(2); iZ <= V.Vy.fCore.ubound(2); iZ += iX) {
                xCoord = mesh.xStaggr(V.Vy.fCore.lbound(0)) - (mesh.xStaggr(V.Vy.fCore.lbound(0) + iX) - mesh.xStaggr(V.Vy.fCore.lbound(0)));
                nseRHS.Vy(V.Vy.fCore.lbound(0) - iX, iY, iZ) = -cos(2.0*M_PI*xCoord/mesh.xLen)*
                                                                sin(2.0*M_PI*mesh.yColloc(iY)/mesh.yLen)*
                                                                cos(2.0*M_PI*mesh.zStaggr(iZ)/mesh.zLen);

                xCoord = mesh.xStaggr(V.Vy.fCore.ubound(0)) + (mesh.xStaggr(V.Vy.fCore.ubound(0)) - mesh.xStaggr(V.Vy.fCore.ubound(0) - iX));
                nseRHS.Vy(V.Vy.fCore.ubound(0) + iX, iY, iZ) = -cos(2.0*M_PI*xCoord/mesh.xLen)*
                                                                sin(2.0*M_PI*mesh.yColloc(iY)/mesh.yLen)*
                                                                cos(2.0*M_PI*mesh.zStaggr(iZ)/mesh.zLen);
            }
        }
    }

    // Y-VELOCITY IN FRONT AND BACK PADS
    for (int iY = 1; iY <= mesh.padWidths(1); iY++) {
        for (int iX = V.Vy.fCore.lbound(0); iX <= V.Vy.fCore.ubound(0); iX += iY) {
            for (int iZ = V.Vy.fCore.lbound(2); iZ <= V.Vy.fCore.ubound(2); iZ += iY) {
                yCoord = mesh.yColloc(V.Vy.fCore.lbound(1)) - (mesh.yColloc(V.Vy.fCore.lbound(1) + iY) - mesh.yColloc(V.Vy.fCore.lbound(1)));
                nseRHS.Vy(iX, V.Vy.fCore.lbound(1) - iY, iZ) = -cos(2.0*M_PI*mesh.xStaggr(iX)/mesh.xLen)*
                                                                sin(2.0*M_PI*yCoord/mesh.yLen)*
                                                                cos(2.0*M_PI*mesh.zStaggr(iZ)/mesh.zLen);

                yCoord = mesh.yColloc(V.Vy.fCore.ubound(1)) + (mesh.yColloc(V.Vy.fCore.ubound(1)) - mesh.yColloc(V.Vy.fCore.ubound(1) - iY));
                nseRHS.Vy(iX, V.Vy.fCore.ubound(1) + iY, iZ) = -cos(2.0*M_PI*mesh.xStaggr(iX)/mesh.xLen)*
                                                                sin(2.0*M_PI*yCoord/mesh.yLen)*
                                                                cos(2.0*M_PI*mesh.zStaggr(iZ)/mesh.zLen);
            }
        }
    }

    // Y-VELOCITY IN BOTTOM AND TOP PADS
    for (int iZ = 1; iZ <= mesh.padWidths(2); iZ++) {
        for (int iX = V.Vy.fCore.lbound(0); iX <= V.Vy.fCore.ubound(0); iX += iZ) {
            for (int iY = V.Vy.fCore.lbound(1); iY <= V.Vy.fCore.ubound(1); iY += iZ) {
                zCoord = mesh.zStaggr(V.Vy.fCore.lbound(2)) - (mesh.zStaggr(V.Vy.fCore.lbound(2) + iZ) - mesh.zStaggr(V.Vy.fCore.lbound(2)));
                nseRHS.Vy(iX, iY, V.Vy.fCore.lbound(2) - iZ) = -cos(2.0*M_PI*mesh.xStaggr(iX)/mesh.xLen)*
                                                                sin(2.0*M_PI*mesh.yColloc(iY)/mesh.yLen)*
                                                                cos(2.0*M_PI*zCoord/mesh.zLen);

                zCoord = mesh.zStaggr(V.Vy.fCore.ubound(2)) + (mesh.zStaggr(V.Vy.fCore.ubound(2)) - mesh.zStaggr(V.Vy.fCore.ubound(2) - iZ));
                nseRHS.Vy(iX, iY, V.Vy.fCore.ubound(2) + iZ) = -cos(2.0*M_PI*mesh.xStaggr(iX)/mesh.xLen)*
                                                                sin(2.0*M_PI*mesh.yColloc(iY)/mesh.yLen)*
                                                                cos(2.0*M_PI*zCoord/mesh.zLen);
            }
        }
    }

    V.imposeBCs();

    V -= nseRHS;

    return std::max(blitz::max(fabs(V.Vx.F)), blitz::max(fabs(V.Vy.F)));
}


scalar_d3::~scalar_d3() { }
