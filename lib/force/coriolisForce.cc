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
/*! \file coriolisForce.cc
 *
 *  \brief Definitions for functions of class coriolisForce
 *  \sa force.h
 *  \author Shashwat Bhattacharya, Roshan Samuel
 *  \date Nov 2019
 *  \copyright New BSD License
 *
 ********************************************************************************************************************************************
 */

#include "force.h"

coriolisForce::coriolisForce(const grid &mesh, vfield &U): force(mesh, U) {
    Fr = 1.0/mesh.inputParams.Ro;
}


void coriolisForce::addForcing(plainvf &Hv) {
    /*
    //ADD THE ROTATING TERM TO THE Vx COMPONENT OF Hv
    V.interTempX = 0.0;
    for (unsigned int i=0; i < V.Vx.VyIntSlices.size(); i++) {
        V.interTempX(V.Vx.fCore) += V.Vy.F(V.Vx.VyIntSlices(i));
    }   
    V.interTempX /= V.Vx.VyIntSlices.size();

    Hv.Vx += Fr*V.interTempX;

    //SUBTRACT THE ROTATING TERM FROM THE Vy COMPONENT of Hv
    V.interTempY = 0.0;
    for (unsigned int i=0; i < V.Vy.VxIntSlices.size(); i++) {
        V.interTempY(V.Vy.fCore) += V.Vx.F(V.Vy.VxIntSlices(i));
    }   
    V.interTempY /= V.Vy.VxIntSlices.size();

    Hv.Vy -= Fr*V.interTempY;
    */
    
    
       //*******Rotation in -x-direction. Edited by Md. Anas *******//

    //(2Omegax uz ey)ADD THE ROTATING TERM TO THE Vy COMPONENT OF Hv
    V.interTempY = 0.0;
    for (unsigned int i=0; i < V.Vy.VzIntSlices.size(); i++) {
        V.interTempY(V.Vy.fCore) += V.Vz.F(V.Vy.VzIntSlices(i));
    }
    V.interTempY /= V.Vy.VzIntSlices.size();

    //Hv.Vy += Fr*V.interTempY;
    Hv.Vy -= Fr*V.interTempY;


    //(-2Omegax uy ez) SUBTRACT THE ROTATING TERM FROM THE Vz COMPONENT of Hv
    V.interTempZ = 0.0;
    for (unsigned int i=0; i < V.Vz.VyIntSlices.size(); i++) {
        V.interTempZ(V.Vz.fCore) += V.Vy.F(V.Vz.VyIntSlices(i));
    }
    V.interTempZ /= V.Vz.VyIntSlices.size();

    //Hv.Vz -= Fr*V.interTempZ;
    Hv.Vz += Fr*V.interTempZ;
    
}
