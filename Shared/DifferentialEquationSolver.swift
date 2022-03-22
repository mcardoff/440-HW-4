//
//  DifferentialEquationSolver.swift
//  Schrodinger-Solver
//
//  Created by Michael Cardiff on 2/18/22.
//

/*
 d
 -- / x(t) \ = / f(x,y,t) \
 dt \ y(t) /   \ g(x,y,t) /
 
 for the Schrodinger equation:
 d(psi)
 ------ = psip
   dt

 d(psip)      m (E - V)
 ------- = -2 --------- psi
    dt          hb^2
 
 linearize derivative:
 
 psi(x+h) - psi(x)
 ----------------- = psip(x) => psi(x+h) = psi(x) + h (psip(x))
          h

 psip(x+h) - psip(x)      m (E - V)                                    /     m (E - V(x)) \
 ------------------- = -2 --------- psi(x) => psip(x+h) = psip(x) + h  | -2 ------------- |
           h                 hb^2                                      \       hb^2       /
 
 */

import Foundation

typealias PotentialFunc = (_:Double) -> Double
typealias PotentialList = (xs: [Double], Vs: [Double])
typealias InitialCondition = (psi: Double, psip: Double)
typealias PhaseSpacePt = (psi: Double, psip: Double)
typealias Iterfunctype = (_: Double, _: Double, _: Double, _: Double, _: Double, _: Double) -> PhaseSpacePt

class SchrodingerSolver: NSObject, ObservableObject {
    
    @Published var totalFuncToPlot : [[plotDataType]] = []
    @Published var energyFunctional : [plotDataType] = []
    @Published var potentialPlot : [plotDataType] = []
    @Published var energyEigenValues : [Double] = []
    @Published var eigenvalueText : String = ""
    
    /// boundaryValProblem
    /// Solve the Boundary Value problem associated with the Schrodinger equation and its potential
    /// - Parameters:
    ///   - a: The Well Width
    ///   - steps: number of points to put in between 0,a
    ///   - Vf: The potential function
    ///   - ic: the initial condition
    ///   - eMin: Min energy to search
    ///   - eMax: Max energy to search
    ///   - eStride: Step size to search for energy eigenvalues
    /// - returns: Normalized wavefunctions
    func boundaryValProblem(a: Double, steps: Int, Vf: PotentialList, ic: InitialCondition,
                            eMin: Double, eMax: Double, eStride: Double) {
        let precision = 1e-9
        // these will have the good energy eigenvalues and the functions as well:
        var goodEnergyPsiCollection : [[Double]] = [], goodEnergyValCollection : [Double] = []
        var energyFunc : [(psi: Double, energy: Double)] = []
        
        // output from rknSolve
//        let lastPointForBC : [(psi: Double, energy: Double)] = rk4Solve(a: a, steps: steps, Vf: Vf, ic: ic,
//                                                                        eMin: eMin, eMax: eMax, eStride: eStride)
        var prevEnergy : Double = 0.0, prevPsi = 0.0
        for energyVal in stride(from: eMin, to: eMax, by: eStride) {
            let thing = rk4SingleEigenVal(a: a, steps: steps, energyVal: energyVal, Vf: Vf, ic: ic)
            let psiVal = thing.lastVal
            
            energyFunc.append((psi: psiVal, energy: energyVal))
            
            // check sign of function output dumbass
            if (prevPsi.sign != psiVal.sign) {
                var checkedPsi = psiVal
                // recalculate the value of the functional with this energy
                var leftVal = prevEnergy, rightVal = energyVal
                var testVal = 0.0
                var possibleAnswer : (totalPsi: [Double], lastVal: Double) = (totalPsi: [], -12.6)

                while(abs(checkedPsi) > precision) {
                    // method of secants
                    testVal = rightVal - checkedPsi * (rightVal - leftVal) / (checkedPsi - prevPsi)
                    possibleAnswer = rk4SingleEigenVal(a: a, steps: steps, energyVal: testVal, Vf: Vf, ic: ic)
                    let possibleZero = possibleAnswer.lastVal
                    
                    leftVal = rightVal
                    rightVal = testVal
                    
                    prevPsi = checkedPsi
                    checkedPsi = possibleZero
                    energyFunc.append((psi: possibleZero, testVal))
                }
                // now the value is in a good range
                if(!possibleAnswer.totalPsi.isEmpty) {
                    goodEnergyPsiCollection.append(possibleAnswer.totalPsi)
                    energyEigenValues.append(testVal)
                }
            }
            
            prevEnergy = energyVal
            prevPsi = psiVal
            energyFunc.append((psi: psiVal, energy: energyVal))
        }
        
        var tempStringList : [String] = []
        for num in energyEigenValues {
            tempStringList.append(String(format: "%0.16f", num))
        }
        
        eigenvalueText = tempStringList.joined(separator: "\n")
        
        // make sure the plotting does not look weird
        energyFunc.sort(by: {
            (x: (psi: Double, energy: Double), y: (psi: Double, energy: Double)) -> Bool in return x.energy < y.energy
        })
        
        goodEnergyPsiCollection = normalizeWaveFuncs(psiCollection: goodEnergyPsiCollection, a: a, steps: steps)
        
        toPlotData(xvals: Vf.xs, yvals: goodEnergyPsiCollection)
        energyEigenValues.append(contentsOf: goodEnergyValCollection)
        fillPotentialPlot(potential: Vf)
        fillEnergyFunc(vals: energyFunc)
    }
    
    /// normalizeWaveFuncs
    /// use the average value theorem to normalize the wavefunctions
    /// - Parameters:
    ///   - psiCollection: functions to normalize
    ///   - a: The Well Width
    ///   - steps: number of points to put in between 0,a
    /// - returns: Normalized wavefunctions
    func normalizeWaveFuncs(psiCollection: [[Double]], a: Double, steps: Int) -> [[Double]] {
        // integrate the function and use thact factor to ensure normalization
        assert(psiCollection.count > 0)
        var newPsiCollection : [[Double]] = []
        for list in psiCollection {
//            assert(list.count == steps)
            
            var sum = 0.0
            for item in list {
                sum += item * item
            }
            // average value theorem: I(a->b) = (b - a) * <f>
            let normVal = (a - 0.0) * sum / Double(steps)
            var newList : [Double] = []
            for item in list {
                newList.append(item / normVal)
            }
            newPsiCollection.append(newList)
        }
        return newPsiCollection
    }
    
    /// rknSingleEigenVal:
    /// Solve the 1D Schrodinger Equation for a single energy Eigenvalue and using any iteration function
    /// - Parameters:
    ///   - a: Box width
    ///   - steps: number of points to put in between 0,a
    ///   - energyVal: energy eigenvalue to solve the schrodinger equation with
    ///   - Vf: Potential to solve for
    ///   - ic: Initial Condition
    ///   - iterfunc: iteration function to use, either Euler's method or RK4
    /// - returns: The wavefunction and its last value for use in solving the boundary value problem
    func rknSingleEigenVal(a: Double, steps: Int, energyVal: Double, Vf: PotentialList, ic: InitialCondition, iterfunc: Iterfunctype) -> (totalPsi: [Double], lastVal: Double) {
        let xs : [Double] = Vf.xs, vs : [Double] = Vf.Vs
        let stepSize = a / Double(steps), h = stepSize
        var curPsi = ic.psi, curPsip = ic.psip
        var psiList : [Double] = []
        for (x,V) in zip(xs, vs) {
//            let V = vs[i], x = xs[i]
            
            let deltatup : PhaseSpacePt = iterfunc(h, curPsi, curPsip, x, energyVal, V)
            
            let nextPsi = curPsi + deltatup.psi,
                nextPsiP = curPsip + deltatup.psip
            
            curPsi  = nextPsi
            curPsip = nextPsiP
            
            psiList.append(curPsi)
        }
        return (totalPsi: psiList, lastVal: curPsi)
    }
    
    /// rk4SingleEigenVal:
    /// Solve the 1D Schrodinger Equation for a single energy Eigenvalue using rk4
    /// - Parameters:
    ///   - a: Box width
    ///   - steps: number of points to put in between 0,a
    ///   - energyVal: energy eigenvalue to solve the schrodinger equation with
    ///   - Vf: Potential to solve for
    ///   - ic: Initial Condition
    /// - returns: The wavefunction and its last value for use in solving the boundary value problem
    func rk4SingleEigenVal(a: Double, steps: Int, energyVal: Double, Vf: PotentialList, ic: InitialCondition) -> (totalPsi: [Double], lastVal: Double) {
        /// rk4
        /// Matching iterfunc for the below function
        func rk4(h: Double, psi: Double, psip: Double, x: Double, energy: Double, V: Double) -> PhaseSpacePt {
            let k0 = h * psiIter(psi: psi, psip: psip, xval: x),
                l0 = h * psiPrimeIter(psi: psi, psip: psip, xval: x, energy: energy, V: V)
            
            let k1 = h * psiIter(psi: psi+(k0/2), psip: psip+(l0/2), xval: x+(h/2)),
                l1 = h * psiPrimeIter(psi: psi+(k0/2), psip: psip+(l0/2), xval: x+(h/2), energy: energy, V: V)
            
            let k2 = h * psiIter(psi: psi+(k1/2), psip: psip+(l1/2), xval: x+(h/2)),
                l2 = h * psiPrimeIter(psi: psi+(k1/2), psip: psip+(l1/2), xval: x+(h/2), energy: energy, V: V)
            
            let k3 = h * psiIter(psi: psi+h, psip: psip+h, xval: x+h),
                l3 = h * psiPrimeIter(psi: psi+h, psip: psip+h, xval: x+h, energy: energy, V: V)
            
            return (psi: (k0 + 2*k1 + 2*k2 + k3)/6, psip: (l0 + 2*l1 + 2*l2 + l3)/6)
        }
        
        return rknSingleEigenVal(a: a, steps: steps, energyVal: energyVal, Vf: Vf, ic: ic, iterfunc: rk4)
    }
    
    /// psiIter
    /// Functions aiding the RK4 solving method, in order to make code more reuseable
    func psiIter(psi: Double, psip: Double, xval: Double) -> Double {
        return psip
    }
    
    /// psiPrimeIter
    /// Functions aiding the RK4 solving method, in order to make code more reuseable
    func psiPrimeIter(psi: Double, psip: Double, xval: Double, energy: Double, V: Double) -> Double {
        return Schrod(x: xval, mass: 1.0, hbar: 1.0, energy: energy, V: V) * psi
    }
    
    /// Schrod
    /// Specific implementation of the schrodinger equation when turned into a system and solved for the derivative term:
    func Schrod(x: Double, mass: Double, hbar: Double, energy: Double, V: Double) -> Double {
        let consts = -2.0 * mass / (hbar * hbar)
        let rest = energy - V
        return consts * rest
    }
    
    /// Convert arrays of doubles to plotDataTypes
    func toPlotData(xvals: [Double], yvals: [[Double]]) {
        totalFuncToPlot = []
        for vals in yvals {
            var tempList : [plotDataType] = []
            for (x,y) in zip(xvals, vals) {
                tempList.append([.X: x, .Y: y])
            }
            totalFuncToPlot.append(tempList)
        }
    }
    
    /// Create the plotDataType with the energy functional values
    func fillEnergyFunc(vals: [(psi: Double, energy: Double)]) {
        for (psi, energy) in vals {
            energyFunctional.append([.X: energy, .Y: psi])
        }
    }
    
    /// Create the plotDataType with the Potential Values
    func fillPotentialPlot(potential: PotentialList) {
        for (x, V) in zip(potential.xs, potential.Vs) {
            potentialPlot.append([.X: x, .Y: V])
        }
    }
    
    /// Clear the datya
    func clearData() {
        self.totalFuncToPlot.removeAll()
        self.energyFunctional.removeAll()
        self.potentialPlot.removeAll()
        self.energyEigenValues.removeAll()
    }
    
    
    // The following is all deprecated code
    func psiPrimeIterFunctionV(psi: Double, psip: Double, xval: Double, energy: Double, V: PotentialFunc) -> Double {
        return SchrodFunctionV(x: xval, mass: 1.0, hbar: 1.0, energy: energy, V: V) * psi
    }
    
    func SchrodFunctionV(x: Double, mass: Double, hbar: Double, energy: Double, V: PotentialFunc) -> Double {
        let consts = -2.0 * mass / (hbar * hbar)
        let rest = energy - V(x)
        return consts * rest
    }
    
    /// rknSolve:
    /// Solve the 1D Schrodinger Equation using any function
    /// - Parameters:
    ///   - a: Box width
    ///   - steps: number of points to put in between 0,a
    ///   - Vf: Potential to solve for
    ///   - ic: Initial Condition
    ///   - iterfunc: determines which method of solution, euler for RK0, rk4 for RK4
    func rknSolve(a: Double, steps: Int, Vf: PotentialList, ic: InitialCondition, iterfunc: Iterfunctype,
                  eMin: Double, eMax: Double, eStride: Double) -> [(psi: Double, energy: Double)] {
        // Important for the iteration
        let xs : [Double] = Vf.xs, vs : [Double] = Vf.Vs
        let stepSize = a / Double(steps), h = stepSize
        var psiCollection : [[Double]] = [], psipCollection : [[Double]] = []
        var curPsi = ic.psi, curPsip = ic.psip
        
        // root finding
        var lastPointForBC : [(psi: Double, energy: Double)] = []
        
        for energyVal in stride(from: eMin, to: eMax, by: eStride) {
//        let energyVal = Double.pi*Double.pi / 8
//            print(energyVal)
            var curPsiList = [ic.psi], curPsipList = [ic.psip]
            for i in 1..<xs.count {
                let V = vs[i], x = xs[i]
                
                let deltatup : PhaseSpacePt = iterfunc(h, curPsi, curPsip, x, energyVal, V)
                
                let nextPsi = curPsi + deltatup.psi,
                    nextPsiP = curPsip + deltatup.psip
                
                curPsi  = nextPsi
                curPsip = nextPsiP
                
                curPsiList.append(curPsi)
                curPsipList.append(curPsip)
            }
            
            psiCollection.append(curPsiList)
            psipCollection.append(curPsipList)
            lastPointForBC.append((psi: curPsi, energy: energyVal))
            
            // reset the state dumbass
            curPsi = ic.psi
            curPsip = ic.psip

        }
        
        return lastPointForBC
    }
    
    /// eulerSolve:
    /// Solve the 1D Schrodinger Equation using euler's method
    /// - Parameters:
    ///   - a: Box width
    ///   - steps: number of points to put in between 0,a
    ///   - Vf: Potential to solve for
    ///   - ic: Initial Condition
    func eulerSolve(a: Double, steps: Int, Vf: PotentialList, ic: InitialCondition,
                    eMin: Double, eMax: Double, eStride: Double) -> [(psi: Double, energy: Double)] {
        /// euler
        /// Matching iterfunc for the parent function
        func euler(h: Double, psi: Double, psip: Double, x: Double, energy: Double, V: Double) -> PhaseSpacePt {
            let k0 = h * psiIter(psi: psi, psip: psip, xval: x),
                l0 = h * psiPrimeIter(psi: psi, psip: psip, xval: x, energy: energy, V: V)
            
            return (psi: k0, psip: l0)
        }
        
        return rknSolve(a: a, steps: steps, Vf: Vf, ic: ic, iterfunc: euler, eMin: eMin, eMax: eMax, eStride: eStride)
    }
    
    /// rk4Solve:
    /// Solve the 1D Schrodinger Equation using rk4
    /// - Parameters:
    ///   - a: Box width
    ///   - steps: number of points to put in between 0,a
    ///   - Vf: Potential to solve for
    ///   - ic: Initial Condition
    func rk4Solve(a: Double, steps: Int, Vf: PotentialList, ic: InitialCondition,
                  eMin: Double, eMax: Double, eStride: Double) -> [(psi: Double, energy: Double)] {
        /// rk4
        /// Matching iterfunc for the above function
        func rk4(h: Double, psi: Double, psip: Double, x: Double, energy: Double, V: Double) -> PhaseSpacePt {
            let k0 = h * psiIter(psi: psi, psip: psip, xval: x),
                l0 = h * psiPrimeIter(psi: psi, psip: psip, xval: x, energy: energy, V: V)
            
            let k1 = h * psiIter(psi: psi+(k0/2), psip: psip+(l0/2), xval: x+(h/2)),
                l1 = h * psiPrimeIter(psi: psi+(k0/2), psip: psip+(l0/2), xval: x+(h/2), energy: energy, V: V)
            
            let k2 = h * psiIter(psi: psi+(k1/2), psip: psip+(l1/2), xval: x+(h/2)),
                l2 = h * psiPrimeIter(psi: psi+(k1/2), psip: psip+(l1/2), xval: x+(h/2), energy: energy, V: V)
            
            let k3 = h * psiIter(psi: psi+h, psip: psip+h, xval: x+h),
                l3 = h * psiPrimeIter(psi: psi+h, psip: psip+h, xval: x+h, energy: energy, V: V)
            
            return (psi: (k0 + 2*k1 + 2*k2 + k3)/6, psip: (l0 + 2*l1 + 2*l2 + l3)/6)
        }
        return rknSolve(a: a, steps: steps, Vf: Vf, ic: ic, iterfunc: rk4, eMin: eMin, eMax: eMax, eStride: eStride)
    }
}
