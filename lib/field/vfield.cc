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
/*! \file vfield.cc
 *
 *  \brief Definitions for functions of class vfield - vector field
 *  \sa vfield.h
 *  \author Roshan Samuel, Ali Asad
 *  \date Nov 2019
 *  \copyright New BSD License
 *
 ********************************************************************************************************************************************
 */

#include "plainsf.h"
#include "plainvf.h"
#include "sfield.h"
#include "vfield.h"
#include <math.h>

/**
 ********************************************************************************************************************************************
 * \brief   Constructor of the vfield class
 *
 *          Three instances of the field class to store the data of the three components are initialized.
 *          The fields are initialized with appropriate grid staggering to place the components on the cell faces.
 *          The name for the vector field as given by the user is also assigned.
 *
 * \param   gridData is a const reference to the global data contained in the grid class
 * \param   fieldName is a string value set by the user to name and identify the vector field
 ********************************************************************************************************************************************
 */
vfield::vfield(const grid &gridData, std::string fieldName):
               gridData(gridData),
               Vx(gridData, "Vx", false, true, true),
               Vy(gridData, "Vy", true, false, true),
               Vz(gridData, "Vz", true, true, false),
               derVx(gridData, Vx), derVy(gridData, Vy), derVz(gridData, Vz)
{
    this->fieldName = fieldName;

    interTempX.resize(Vx.fSize);
    interTempX.reindexSelf(Vx.flBound);

    derivTempX.resize(Vx.fSize);
    derivTempX.reindexSelf(Vx.flBound);

#ifndef PLANAR
    interTempY.resize(Vy.fSize);
    interTempY.reindexSelf(Vy.flBound);

    derivTempY.resize(Vy.fSize);
    derivTempY.reindexSelf(Vy.flBound);
#endif

    interTempZ.resize(Vz.fSize);
    interTempZ.reindexSelf(Vz.flBound);

    derivTempZ.resize(Vz.fSize);
    derivTempZ.reindexSelf(Vz.flBound);
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to compute the diffusion term
 *
 *          It is assumed that the vector field is specified at face-centers, as required by the \ref sfield#computeNLin
 *          "computeNLin" function of sfield.
 *
 * \param   H is a pointer to a vector field (vfield) to which the output of the function is to be written
 ********************************************************************************************************************************************
 */
void vfield::computeDiff(plainvf &H) {
    derivTempX = 0.0;
    derVx.calcDerivative2xx(derivTempX);
    H.Vx(Vx.fCore) += derivTempX(Vx.fCore);

#ifndef PLANAR
    derivTempX = 0.0;
    derVx.calcDerivative2yy(derivTempX);
    H.Vx(Vx.fCore) += derivTempX(Vx.fCore);
#endif

    derivTempX = 0.0;
    derVx.calcDerivative2zz(derivTempX);
    H.Vx(Vx.fCore) += derivTempX(Vx.fCore);

#ifndef PLANAR
    derivTempY = 0.0;
    derVy.calcDerivative2xx(derivTempY);
    H.Vy(Vy.fCore) += derivTempY(Vy.fCore);

    derivTempY = 0.0;
    derVy.calcDerivative2yy(derivTempY);
    H.Vy(Vy.fCore) += derivTempY(Vy.fCore);

    derivTempY = 0.0;
    derVy.calcDerivative2zz(derivTempY);
    H.Vy(Vy.fCore) += derivTempY(Vy.fCore);
#endif

    derivTempZ = 0.0;
    derVz.calcDerivative2xx(derivTempZ);
    H.Vz(Vz.fCore) += derivTempZ(Vz.fCore);

#ifndef PLANAR
    derivTempZ = 0.0;
    derVz.calcDerivative2yy(derivTempZ);
    H.Vz(Vz.fCore) += derivTempZ(Vz.fCore);
#endif

    derivTempZ = 0.0;
    derVz.calcDerivative2zz(derivTempZ);
    H.Vz(Vz.fCore) += derivTempZ(Vz.fCore);
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to compute the convective derivative of the vector field
 *
 *          The function computes for the operator \f$ (\mathbf{u}.\nabla)\mathbf{v} \f$ for the vector field, \f$\mathbf{v}\f$.
 *          To do so, the function needs the vector field (vfield) of velocity, \f$\mathbf{u}\f$.
 *          This value is used in three separate calls to the \ref sfield#computeNLin "computeNLin" function of sfield
 *          to compute the derivatives for the three components of the vector field.
 *          It is assumed that the velocity is specified at face-centers, as required by the \ref sfield#computeNLin
 *          "computeNLin" function of sfield.
 *
 * \param   V is a const reference to a vector field (vfield) that specifies the convection velocity at each point
 * \param   H is a pointer to a vector field (vfield) to which the output of the function is to be written
 ********************************************************************************************************************************************
 */
void vfield::computeNLin(const vfield &V, plainvf &H) {
    // Compute non-linear term for the Vx component
    interTempX = 0.0;
    for (unsigned int i=0; i < Vx.VxIntSlices.size(); i++) {
        interTempX(Vx.fCore) += V.Vx.F(Vx.VxIntSlices(i));
    }

    derivTempX = 0.0;
    derVx.calcDerivative1_x(derivTempX);
    H.Vx(Vx.fCore) -= interTempX(Vx.fCore)*derivTempX(Vx.fCore)/Vx.VxIntSlices.size();

#ifndef PLANAR
    interTempX = 0.0;
    for (unsigned int i=0; i < Vx.VyIntSlices.size(); i++) {
        interTempX(Vx.fCore) += V.Vy.F(Vx.VyIntSlices(i));
    }

    derivTempX = 0.0;
    derVx.calcDerivative1_y(derivTempX);
    H.Vx(Vx.fCore) -= interTempX(Vx.fCore)*derivTempX(Vx.fCore)/Vx.VyIntSlices.size();
#endif

    interTempX = 0.0;
    for (unsigned int i=0; i < Vx.VzIntSlices.size(); i++) {
        interTempX(Vx.fCore) += V.Vz.F(Vx.VzIntSlices(i));
    }

    derivTempX = 0.0;    
    derVx.calcDerivative1_z(derivTempX);
    H.Vx(Vx.fCore) -= interTempX(Vx.fCore)*derivTempX(Vx.fCore)/Vx.VzIntSlices.size();

// Compute non-linear term for the Vy component
#ifndef PLANAR
    interTempY = 0.0;
    for (unsigned int i=0; i < Vy.VxIntSlices.size(); i++) {
        interTempY(Vy.fCore) += V.Vx.F(Vy.VxIntSlices(i));
    }

    derivTempY = 0.0;
    derVy.calcDerivative1_x(derivTempY);
    H.Vy(Vy.fCore) -= interTempY(Vy.fCore)*derivTempY(Vy.fCore)/Vy.VxIntSlices.size();

    interTempY = 0.0;
    for (unsigned int i=0; i < Vy.VyIntSlices.size(); i++) {
        interTempY(Vy.fCore) += V.Vy.F(Vy.VyIntSlices(i));
    }

    derivTempY = 0.0;
    derVy.calcDerivative1_y(derivTempY);
    H.Vy(Vy.fCore) -= interTempY(Vy.fCore)*derivTempY(Vy.fCore)/Vy.VyIntSlices.size();

    interTempY = 0.0;
    for (unsigned int i=0; i < Vy.VzIntSlices.size(); i++) {
        interTempY(Vy.fCore) += V.Vz.F(Vy.VzIntSlices(i));
    }

    derivTempY = 0.0;
    derVy.calcDerivative1_z(derivTempY);
    H.Vy(Vy.fCore) -= interTempY(Vy.fCore)*derivTempY(Vy.fCore)/Vy.VzIntSlices.size();
#endif

    // Compute non-linear term for the Vz component
    interTempZ = 0.0;
    for (unsigned int i=0; i < Vz.VxIntSlices.size(); i++) {
        interTempZ(Vz.fCore) += V.Vx.F(Vz.VxIntSlices(i));
    }

    derivTempZ = 0.0;
    derVz.calcDerivative1_x(derivTempZ);
    H.Vz(Vz.fCore) -= interTempZ(Vz.fCore)*derivTempZ(Vz.fCore)/Vz.VxIntSlices.size();

#ifndef PLANAR
    interTempZ = 0.0;
    for (unsigned int i=0; i < Vz.VyIntSlices.size(); i++) {
        interTempZ(Vz.fCore) += V.Vy.F(Vz.VyIntSlices(i));
    }

    derivTempZ = 0.0;
    derVz.calcDerivative1_y(derivTempZ);
    H.Vz(Vz.fCore) -= interTempZ(Vz.fCore)*derivTempZ(Vz.fCore)/Vz.VyIntSlices.size();
#endif

    interTempZ = 0.0;
    for (unsigned int i=0; i < Vz.VzIntSlices.size(); i++) {
        interTempZ(Vz.fCore) += V.Vz.F(Vz.VzIntSlices(i));
    }

    derivTempZ = 0.0;
    derVz.calcDerivative1_z(derivTempZ);
    H.Vz(Vz.fCore) -= interTempZ(Vz.fCore)*derivTempZ(Vz.fCore)/Vz.VzIntSlices.size();
}

/**
 ********************************************************************************************************************************************
 * \brief   Operator is used to calculate time step #dt_out for time integration using CFL Condition with desired Courant No
 *           
 * \param   dt_out is a reference to a real variable into which the calculated value of time-step has to be written
 *********************************************************************************************************************************************
 */
void vfield::computeTStp(real &dt_out) {
    real Umax, Vmax, Wmax;
    real delx, dely, delz; 

    delx = gridData.dXi;
    Umax = Vx.fieldMax();
#ifdef PLANAR
    dely = 0.0;
    Vmax = 1.0;
#else
    dely = gridData.dEt;
    Vmax = Vy.fieldMax();
#endif
    delz = gridData.dZt;
    Wmax = Vz.fieldMax();

    dt_out = real(gridData.inputParams.courantNumber*std::min(std::min(delx/Umax, dely/Vmax), std::min(dely/Vmax, delz/Wmax)));

}

/**
 ********************************************************************************************************************************************
 * \brief   Operator to compute the divergence of the vector field
 *
 *          The operator computes the divergence of a face-centered staggered vector field, and stores it into a cell centered
 *          scalar field as defined by the tensor operation:
 *          \f$ \nabla . \mathbf{v} = \frac{\partial \mathbf{v}}{\partial x} +
 *                                    \frac{\partial \mathbf{v}}{\partial y} +
 *                                    \frac{\partial \mathbf{v}}{\partial z} \f$.
 *
 * \param   divV is a pointer to a scalar field (sfield) into which the computed divergence must be written.
 ********************************************************************************************************************************************
 */
void vfield::divergence(plainsf &divV, const sfield &P) {
    blitz::firstIndex i;
    blitz::secondIndex j;
    blitz::thirdIndex k;

    blitz::Range xStag, yStag, zStag;

    xStag = blitz::Range(P.F.fCore.lbound(0), P.F.fCore.ubound(0));
    yStag = blitz::Range(P.F.fCore.lbound(1), P.F.fCore.ubound(1));
    zStag = blitz::Range(P.F.fCore.lbound(2), P.F.fCore.ubound(2));

    divV = 0.0;

#ifdef PLANAR
    divV.F(P.F.fCore) = gridData.xi_xStaggr(xStag)(i)*(Vx.F(P.F.fCore) - Vx.F(P.F.fCLft))/gridData.dXi + 
                        gridData.zt_zStaggr(zStag)(k)*(Vz.F(P.F.fCore) - Vz.F(P.F.fCBot))/gridData.dZt;
#else
    divV.F(P.F.fCore) = gridData.xi_xStaggr(xStag)(i)*(Vx.F(P.F.fCore) - Vx.F(P.F.fCLft))/gridData.dXi + 
                        gridData.et_yStaggr(yStag)(j)*(Vy.F(P.F.fCore) - Vy.F(P.F.fCFrt))/gridData.dEt + 
                        gridData.zt_zStaggr(zStag)(k)*(Vz.F(P.F.fCore) - Vz.F(P.F.fCBot))/gridData.dZt;
#endif
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to synchronise data across all processors when performing parallel computations
 *
 *          Each of the individual scalar field components have their own subroutine, \ref sfield#syncData "syncData" to send and
 *          receive data across its MPI decomposed sub-domains.
 *          This function calls the \ref sfield#syncData "syncData" function of its components to update the sub-domain boundary pads.
 ********************************************************************************************************************************************
 */
void vfield::syncData() {
    Vx.syncData();
    Vy.syncData();
    Vz.syncData();
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to impose the boundary conditions for the X-component of the vector field
 *
 *          The function first calls the \ref field#syncData "syncData" function of the Vx field to update the sub-domain pads.
 *          Then the boundary conditions are applied at the full domain boundaries by calling the imposeBC()
 *          of each boundary class object assigned to each wall.
 *          The order of imposing boundary conditions is - left, right, front, back, bottom and top boundaries.
 *          The corner values are not being imposed specifically and is thus dependent on the above order.
 *
 ********************************************************************************************************************************************
 */
void vfield::imposeVxBC() {
    Vx.syncData();

    if (not gridData.inputParams.xPer) {
        uLft->imposeBC();
        uRgt->imposeBC();
    }
#ifndef PLANAR
    if (not gridData.inputParams.yPer) {
        uFrn->imposeBC();
        uBak->imposeBC();
    }
#endif
    uTop->imposeBC();
    uBot->imposeBC();
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to impose the boundary conditions for the Y-component of the vector field
 *
 *          The function first calls the \ref field#syncData "syncData" function of the Vy field to update the sub-domain pads.
 *          Then the boundary conditions are applied at the full domain boundaries by calling the imposeBC()
 *          of each boundary class object assigned to each wall.
 *          The order of imposing boundary conditions is - left, right, front, back, bottom and top boundaries.
 *          The corner values are not being imposed specifically and is thus dependent on the above order.
 *
 ********************************************************************************************************************************************
 */
void vfield::imposeVyBC() {
    Vy.syncData();

    if (not gridData.inputParams.xPer) {
        vLft->imposeBC();
        vRgt->imposeBC();
    }
#ifndef PLANAR
    if (not gridData.inputParams.yPer) {
        vFrn->imposeBC();
        vBak->imposeBC();
    }
#endif
    vTop->imposeBC();
    vBot->imposeBC();
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to impose the boundary conditions for the Z-component of the vector field
 *
 *          The function first calls the \ref field#syncData "syncData" function of the Vz field to update the sub-domain pads.
 *          Then the boundary conditions are applied at the full domain boundaries by calling the imposeBC()
 *          of each boundary class object assigned to each wall.
 *          The order of imposing boundary conditions is - left, right, front, back, bottom and top boundaries.
 *          The corner values are not being imposed specifically and is thus dependent on the above order.
 *
 ********************************************************************************************************************************************
 */
void vfield::imposeVzBC() {
    Vz.syncData();

    if (not gridData.inputParams.xPer) {
        wLft->imposeBC();
        wRgt->imposeBC();
    }
#ifndef PLANAR
    if (not gridData.inputParams.yPer) {
        wFrn->imposeBC();
        wBak->imposeBC();
    }
#endif
    wTop->imposeBC();
    wBot->imposeBC();
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to impose the boundary conditions for all the components of the vfield
 *
 *          The function merely calls \ref vfield#imposeVxBC "imposeVxBC", \ref vfield#imposeVyBC "imposeVyBC"
 *          and \ref vfield#imposeVzBC "imposeVzBC" functions together.
 *
 ********************************************************************************************************************************************
 */
void vfield::imposeBCs() {
    imposeVxBC();
#ifndef PLANAR
    imposeVyBC();
#endif
    imposeVzBC();
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to add a given plain vector field
 *
 *          The unary operator += adds a given plain vector field to the entire field stored as vfield and returns
 *          a pointer to itself.
 *
 * \param   a is a reference to a plainvf to be deducted from the member fields
 *
 * \return  A pointer to itself is returned by the vector field class to which the operator belongs
 ********************************************************************************************************************************************
 */
vfield& vfield::operator += (plainvf &a) {
    Vx.F += a.Vx;
    Vy.F += a.Vy;
    Vz.F += a.Vz;

    return *this;
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to subtract a given plain vector field
 *
 *          The unary operator -= subtracts a given plain vector field from the entire field stored as vfield and returns
 *          a pointer to itself.
 *
 * \param   a is a reference to a plainvf to be deducted from the member fields
 *
 * \return  A pointer to itself is returned by the vector field class to which the operator belongs
 ********************************************************************************************************************************************
 */
vfield& vfield::operator -= (plainvf &a) {
    Vx.F -= a.Vx;
    Vy.F -= a.Vy;
    Vz.F -= a.Vz;

    return *this;
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to add a given vector field
 *
 *          The unary operator += adds a given vector field to the entire field stored as vfield and returns
 *          a pointer to itself.
 *
 * \param   a is a reference to another vfield to be deducted from the member fields
 *
 * \return  A pointer to itself is returned by the vector field class to which the operator belongs
 ********************************************************************************************************************************************
 */
vfield& vfield::operator += (vfield &a) {
    Vx.F += a.Vx.F;
    Vy.F += a.Vy.F;
    Vz.F += a.Vz.F;

    return *this;
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to subtract a given vector field
 *
 *          The unary operator -= subtracts a given vector field from the entire field stored as vfield and returns
 *          a pointer to itself.
 *
 * \param   a is a reference to another vfield to be deducted from the member fields
 *
 * \return  A pointer to itself is returned by the vector field class to which the operator belongs
 ********************************************************************************************************************************************
 */
vfield& vfield::operator -= (vfield &a) {
    Vx.F -= a.Vx.F;
    Vy.F -= a.Vy.F;
    Vz.F -= a.Vz.F;

    return *this;
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to multiply a scalar value to the vector field
 *
 *          The unary operator *= multiplies a real value to all the fields (Vx, Vy and Vz) stored in vfield and returns
 *          a pointer to itself.
 *
 * \param   a is a real number to be multiplied to the vector field
 *
 * \return  A pointer to itself is returned by the vector field class to which the operator belongs
 ********************************************************************************************************************************************
 */
vfield& vfield::operator *= (real a) {
    Vx.F *= a;
    Vy.F *= a;
    Vz.F *= a;

    return *this;
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to assign a plain vector field to the vector field
 *
 *          The operator = assigns all the three scalar sub-fields of a vfield to all the corresponding fields (Vx, Vy and Vz).
 *
 * \param   a is a plainvf to be assigned to the vector field
 ********************************************************************************************************************************************
 */
void vfield::operator = (plainvf &a) {
    Vx.F = a.Vx;
    Vy.F = a.Vy;
    Vz.F = a.Vz;
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to assign another vector field to the vector field
 *
 *          The operator = assigns all the three scalar sub-fields of a vfield to all the corresponding fields (Vx, Vy and Vz).
 *
 * \param   a is a vfield to be assigned to the vector field
 ********************************************************************************************************************************************
 */
void vfield::operator = (vfield &a) {
    Vx.F = a.Vx.F;
    Vy.F = a.Vy.F;
    Vz.F = a.Vz.F;
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to assign a scalar value to the vector field
 *
 *          The operator = assigns a real value to all the fields (Vx, Vy and Vz) stored in vfield.
 *
 * \param   a is a real number to be assigned to the vector field
 ********************************************************************************************************************************************
 */
void vfield::operator = (real a) {
    Vx.F = a;
    Vy.F = a;
    Vz.F = a;
}
