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
    
    @Published var goodFuncToPlot : [plotDataType] = []
    @Published var totalFuncToPlot : [[plotDataType]] = [[[.X:0.0, .Y:0.0]]]
    @Published var energyFunctional : [plotDataType] = []
    @Published var potentialPlot : [plotDataType] = []
    var energyEigenValues : [Double] = []
    
    func boundaryValProblem(a: Double, steps: Int, Vf: PotentialList, ic: InitialCondition) {
        let precision = 1e-3
        // these will have the good energy eigenvalues and the functions as well:
        var goodEnergyPsiCollection : [[Double]] = [], goodEnergyValCollection : [Double] = []
        var energyFunc : [(psi: Double, energy: Double)] = []
        
        // output from rknSolve
        let lastPointForBC : [(psi: Double, energy: Double)] = rk4Solve(a: a, steps: steps, Vf: Vf, ic: ic, eMin: 1.0, eMax: 32.0, eStride: 0.5)
        
        
        var prevEnergy : Double = 0.0, prevPsi = 0.0
        for tup in lastPointForBC {
            let energyVal = tup.energy
            let psiVal = tup.psi
            
            energyFunc.append((psi: psiVal, energy: energyVal))
            
            // check sign of fucntion output dumbass
            if (prevPsi.sign != psiVal.sign) {
                // sign change detected!
                // find energy value between energyVal and prevEnergy
                var checkedPsi = psiVal
                // recalculate the value of the functional with this energy
                var leftVal = prevEnergy, rightVal = energyVal
                var testVal = 0.0
                var possibleAnswer : (totalPsi: [Double], lastVal: Double) = (totalPsi: [], -12.6)

                while(abs(checkedPsi - prevPsi) > precision) {
                    testVal = rightVal - checkedPsi * (rightVal - leftVal) / (checkedPsi - prevPsi)
//                    midEnergyVal = (rightVal + leftVal) / 2.0
                    possibleAnswer = rknSingleEigenVal(a: a, steps: steps, energyVal: testVal, Vf: Vf, ic: ic, iterfunc: rk4)
                    let possibleZero = possibleAnswer.lastVal
                    
                    // check sign change between possibleZero and checkedVal
                    
//                    if (checkedPsi * possibleZero > 0) { // positive, same sign, zero not in this interval
//                        leftVal = midEnergyVal
//                    } else if (checkedPsi * possibleZero < 0) { // negative, sign change in this interval
//                        rightVal = midEnergyVal
//                    } else if (checkedPsi * possibleZero == 0) {
//
//                    } else { // something went wrong
//                        exit(1001)
//                    }
                    
                    leftVal = rightVal
                    rightVal = testVal
                    
                    prevPsi = checkedPsi
                    checkedPsi = possibleZero
                    energyFunc.append((psi: possibleZero, testVal))
                }
                // now the value is in a good range
                if(!possibleAnswer.totalPsi.isEmpty) {
                    goodEnergyPsiCollection.append(possibleAnswer.totalPsi)
                }
            }
            
            prevEnergy = energyVal
            prevPsi = psiVal
            energyFunc.append((psi: psiVal, energy: energyVal))
        }
        if(goodEnergyPsiCollection.isEmpty) {
            print("EMPTY")
        }
        
        energyFunc.sort(by: secondItem)
        
        toPlotData(xvals: Vf.xs, yvals: goodEnergyPsiCollection)
        energyEigenValues.append(contentsOf: goodEnergyValCollection)
        fillPotentialPlot(potential: Vf)
        fillEnergyFunc(vals: energyFunc)
    }
    
    func secondItem(x: (psi: Double, energy: Double), y: (psi: Double, energy: Double)) -> Bool {
        return x.energy < y.energy
    }
    
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
    
    /// rknSolve:
    /// Solve the 1D Schrodinger Equation using any function
    /// - Parameters:
    ///   - a: Box width
    ///   - steps: number of points to put in between 0,a
    ///   - Vf: Potential to solve for
    ///   - ic: Initial Condition
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
        return rknSolve(a: a, steps: steps, Vf: Vf, ic: ic, iterfunc: rk4, eMin: eMin, eMax: eMax, eStride: eStride)
    }
    
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
    
    /// Helper functions to turn data values in RK4/Euler to plottable things
    func toPlotData(xvals: [Double], yvals: [Double]) {
        for i in 0..<xvals.count {
            let x = xvals[i], y = yvals[i]
            goodFuncToPlot.append([.X: x, .Y: y])
        }
    }
    
    func toPlotData(xvals: [Double], yvals: [[Double]]) {
        totalFuncToPlot = []
        for vals in yvals {
            var tempList : [plotDataType] = []
            for (x,y) in zip(xvals, vals) {
                tempList.append([.X: x, .Y: y])
                goodFuncToPlot.append([.X: x, .Y: y])
            }
            totalFuncToPlot.append(tempList)
        }
    }
    
    func fillEnergyFunc(vals: [(psi: Double, energy: Double)]) {
        for (psi, energy) in vals {
            energyFunctional.append([.X: energy, .Y: psi])
        }
    }
    
    func fillPotentialPlot(potential: PotentialList) {
        for (x, V) in zip(potential.xs, potential.Vs) {
            potentialPlot.append([.X: x, .Y: V])
        }
    }
    
    
    // deprecated versions of the iteration functions for when V was a function not list
    func psiPrimeIterFunctionV(psi: Double, psip: Double, xval: Double, energy: Double, V: PotentialFunc) -> Double {
        return SchrodFunctionV(x: xval, mass: 1.0, hbar: 1.0, energy: energy, V: V) * psi
    }
    
    func SchrodFunctionV(x: Double, mass: Double, hbar: Double, energy: Double, V: PotentialFunc) -> Double {
        let consts = -2.0 * mass / (hbar * hbar)
        let rest = energy - V(x)
        return consts * rest
    }
}
