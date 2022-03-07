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
        return [.square, .linear, .quadratic, .centeredquadratic, .squarebarrier, .trianglebarrier]
    }
    case square
    case linear
    case quadratic
    case centeredquadratic
    case squarebarrier
    case trianglebarrier
    
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
        }
    }
}

func getPotential(xMin: Double, xMax: Double, steps: Int, choice: PotentialType, amplitude: Double) -> PotentialList {
    print(choice)

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
//    default:
//        return squareWell(xMin: xMin, xMax: xMax, steps: steps, height: 0.0)
    }
}

func generalWell(xMin: Double, xMax: Double, steps: Int, f: PotentialFunc) -> PotentialList {
    var x: [Double] = [xMin], V : [Double] = [MXVAL] // start with really high value, should be infinite
    
    let stepSize = (xMax - xMin) / Double(steps)
    for xVal in stride(from: xMin+stepSize, to: xMax-stepSize, by: stepSize) {
        x.append(xVal)
        V.append(f(xVal))
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

func squareBarrier(xMin: Double, xMax: Double, steps: Int, amplitude: Double) -> PotentialList {
    func squareBarrier(x: Double) -> Double {
        if (x > 0.4*xMax && (x < 0.6*xMax)) {
            return amplitude
        } else {
            return 0.0
        }
    }
    
    return generalWell(xMin: xMin, xMax: xMax, steps: steps, f: squareBarrier)
}

func triangleBarrier(xMin: Double, xMax: Double, steps: Int, amplitude: Double) -> PotentialList {
    func triangleBarrier(x: Double) -> Double {
        if (x > 0.4*xMax && (x < xMax / 2.0)) {
            return amplitude * x
        } else if (x > xMax / 2.0 && x < 0.6*xMax) {
            return -amplitude * x + (amplitude * xMax / 2)
        } else {
            return 0.0
        }
    }
    
    return generalWell(xMin: xMin, xMax: xMax, steps: steps, f: triangleBarrier)
}
