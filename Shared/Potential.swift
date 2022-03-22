//
//  Potential.swift
//  Schrodinger-Solver
//
//  Created by Michael Cardiff on 2/18/22.
//

import Foundation

let MXVAL = 100000.0

typealias CoordTuple = (x: Double, y: Double)

enum PotentialType: CaseIterable, Identifiable {
    static var allCases : [PotentialType] {
        return [.square, .linear, .quadratic, .centeredquadratic, .squarebarrier, .trianglebarrier, .squarepluslinear, .coupledQuadratic, .coupledSquarePlusField, .kronigPenney]
    }
    case square
    case linear
    case quadratic
    case centeredquadratic
    case squarebarrier
    case squarepluslinear
    case trianglebarrier
    case coupledQuadratic
    case coupledSquarePlusField
    case kronigPenney
    
    var id: Self { self }
    
    func toString() -> String {
        switch self {
            
        case .square:
            return "Square Well"
        case .linear:
            return "Linear Well"
        case .quadratic:
            return "Quadratic Well"
        case .centeredquadratic:
            return "Centered Quadratic Well"
        case .squarebarrier:
            return "Square Barrier"
        case .trianglebarrier:
            return "Triangle Barrier"
        case .coupledQuadratic:
            return "Coupled Quadratic"
        case .squarepluslinear:
            return "Square+Linear"
        case .coupledSquarePlusField:
            return "Coupled Square+Field"
        case .kronigPenney:
            return "Kronig Penney"
        }
    }
}

func getPotential(xMin: Double, xMax: Double, steps: Int, choice: PotentialType, amplitude: Double) -> PotentialList {
//    print(choice)

    switch(choice) {
    case .square:
        return squareWell(xMin: xMin, xMax: xMax, steps: steps, height: amplitude)
    case .linear:
        return linearWell(xMin: xMin, xMax: xMax, steps: steps, slope: amplitude)
    case .quadratic:
        return quadraticWell(xMin: xMin, xMax: xMax, steps: steps, amplitude: amplitude)
    case .centeredquadratic:
        return centeredQuadraticWell(xMin: xMin, xMax: xMax, steps: steps, amplitude: amplitude)
    case .squarebarrier:
        return squareBarrier(xMin: xMin, xMax: xMax, steps: steps, amplitude: amplitude)
    case .trianglebarrier:
        return triangleBarrier(xMin: xMin, xMax: xMax, steps: steps, amplitude: amplitude)
    case .coupledQuadratic:
        return coupledQuadratic(xMin: xMin, xMax: xMax, steps: steps, amplitude: amplitude)
    case .squarepluslinear:
        return squarePlusLinear(xMin: xMin, xMax: xMax, steps: steps, amplitude: amplitude)
    case .coupledSquarePlusField:
        return coupledSquarePlusField(xMin: xMin, xMax: xMax, steps: steps, amplitude: amplitude)
    case .kronigPenney:
        return kronigPenney(xMin: xMin, xMax: xMax, steps: steps, amplitude: amplitude)
    }
}

func generalWell(xMin: Double, xMax: Double, steps: Int, f: PotentialFunc) -> PotentialList {
    var x: [Double] = [xMin], V : [Double] = [MXVAL] // start with really high value, should be infinite
    
    let stepSize = (xMax - xMin) / Double(steps)
    for xVal in stride(from: xMin+stepSize, to: xMax-stepSize, by: stepSize) {
        x.append(xVal)
        V.append(f(xVal))
//        print(f(xVal))
    }
    
    x.append(xMax)
    V.append(MXVAL)
    return (xs: x, Vs: V)
}

func squareWell(xMin: Double, xMax: Double, steps: Int, height: Double) -> PotentialList {
    return generalWell(xMin: xMin, xMax: xMax, steps: steps, f: {(_:Double) -> Double in return height})
}

func linearWell(xMin: Double, xMax: Double, steps: Int, slope: Double) -> PotentialList {
    return generalWell(xMin: xMin, xMax: xMax, steps: steps, f: {(x:Double) -> Double in return slope * x})
}

func quadraticWell(xMin: Double, xMax: Double, steps: Int, amplitude: Double) -> PotentialList {
    return generalWell(xMin: xMin, xMax: xMax, steps: steps, f: {(x:Double) -> Double in return amplitude * x * x})
}

func squarePlusLinear(xMin: Double, xMax: Double, steps: Int, amplitude: Double) -> PotentialList {
    var x: [Double] = [xMin], V : [Double] = [MXVAL] // start with really high value, should be infinite
    let xStep = (xMax - xMin) / Double(steps)
    
    for i in stride(from: xMin+xStep, to: (xMax+xMin)/2.0, by: xStep) {
        x.append(i)
        V.append(0.0)
    }
    
    for i in stride(from: (xMin+xMax)/2.0, through: xMax-xStep, by: xStep) {
        x.append(i)
        V.append(((i-(xMin+xMax)/2.0)*4.0*0.1))
    }
    
    x.append(xMax)
    V.append(MXVAL)
    
    return (xs: x, Vs: V)
                
}

func centeredQuadraticWell(xMin: Double, xMax: Double, steps: Int, amplitude: Double) -> PotentialList {
    // centered at avg
    func centered(x: Double) -> Double {
        let midpt = (xMin + xMax) / 2
        if (x > 0.5*midpt && (x < 1.5*midpt)) {
            let adjst = (x - midpt)
            return amplitude * adjst * adjst
        } else {
            return 0.0
        }
    }
    
    return generalWell(xMin: xMin, xMax: xMax, steps: steps, f: centered)
}

func coupledQuadratic(xMin: Double, xMax: Double, steps: Int, amplitude: Double) -> PotentialList {
    var x: [Double] = [xMin], V : [Double] = [MXVAL] // start with really high value, should be infinite
    let xStep = (xMax - xMin) / Double(steps)
    
    for i in stride(from: xMin+xStep, to: xMin + (xMax-xMin)*0.5, by: xStep) {
        x.append(i)
        V.append(amplitude*(pow((i-(xMin+(xMax-xMin)/4.0)), 2.0)))
    }
    
    for i in stride(from: xMin + (xMax-xMin)*0.5, through: xMax-xStep, by: xStep) {
        x.append(i)
        V.append(amplitude*(pow((i-(xMax-(xMax-xMin)/4.0)), 2.0)))
    }
    
    x.append(xMax)
    V.append(MXVAL)
    
    return (xs: x, Vs: V)
}

func squareBarrier(xMin: Double, xMax: Double, steps: Int, amplitude: Double) -> PotentialList {
    func squareBarrier(x: Double) -> Double {
        if (x > 0.4*(xMax-xMin) && (x < 0.6*(xMax-xMin))) {
            return amplitude
        } else {
            return 0.0
        }
    }
    
    return generalWell(xMin: xMin, xMax: xMax, steps: steps, f: squareBarrier)
}

func triangleBarrier(xMin: Double, xMax: Double, steps: Int, amplitude: Double) -> PotentialList {
    func triangleBarrier(x: Double) -> Double {
        if ((0.4*xMax < x) && (x < xMax / 2.0)) {
            return amplitude * (x - 0.4*xMax)
        } else if (x >= xMax / 2.0 && x <= 0.6*xMax) {
            return -amplitude * (x - 0.6*xMax)
        } else {
            return 0.0
        }
    }
    
    return generalWell(xMin: xMin, xMax: xMax, steps: steps, f: triangleBarrier)
}

func coupledSquarePlusField(xMin: Double, xMax: Double, steps: Int, amplitude: Double) -> PotentialList {
    var x: [Double] = [xMin], V : [Double] = [MXVAL] // start with really high value, should be infinite
    let xStep = (xMax - xMin) / Double(steps)
    
    for i in stride(from: xMin+xStep, to: xMin + (xMax-xMin)*0.4, by: xStep) {
        
        x.append(i)
        V.append(0.0)
        
    }
    
    for i in stride(from: xMin + (xMax-xMin)*0.4, to: xMin + (xMax-xMin)*0.6, by: xStep) {
        
        x.append(i)
        V.append(amplitude)
        
    }
    
    for i in stride(from: xMin + (xMax-xMin)*0.6, to: xMax, by: xStep) {
        
        x.append(i)
        V.append(0.0)
        
    }
    
    x.append(xMax)
    V.append(MXVAL)
    
    return (xs: x, Vs: V)
}

func kronigPenney(xMin: Double, xMax: Double, steps: Int, amplitude: Double) -> PotentialList {
    var x: [Double] = [xMin], V : [Double] = [MXVAL] // start with really high value, should be infinite
    let xStep = (xMax - xMin) / Double(steps)
    let numberOfBarriers = 2.0
    let barrierPotential = amplitude
    let latticeSpacing = xMax/numberOfBarriers
    let barrierWidth = 1.0/6.0*latticeSpacing
    var currentBarrierPosition = 0.0
    var barrierNumber = 1
    var inBarrier = false
    
    for i in stride(from: xMin+xStep, through: xMax-xStep, by: xStep) {
        
        currentBarrierPosition = -latticeSpacing/2.0 + Double(barrierNumber)*latticeSpacing
        
        if( (abs(i-currentBarrierPosition)) < (barrierWidth/2.0)) {
            inBarrier = true
            x.append(i)
            V.append(barrierPotential)
        }
        else {
            
            if (inBarrier){
                
                inBarrier = false
                barrierNumber += 1
                
            }
            
            x.append(i)
            V.append(0.0)
            
        }
    }
    
    x.append(xMax)
    V.append(MXVAL)
    
    return (xs: x, Vs: V)
}
