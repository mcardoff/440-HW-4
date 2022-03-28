//
//  Tests_macOS.swift
//  Tests macOS
//
//  Created by Michael Cardiff on 2/18/22.
//

import XCTest

class Tests_macOS: XCTestCase {

    override func setUpWithError() throws {
        // Put setup code here. This method is called before the invocation of each test method in the class.

        // In UI tests it is usually best to stop immediately when a failure occurs.
        continueAfterFailure = false

        // In UI tests it’s important to set the initial state - such as interface orientation - required for your tests before they run. The setUp method is a good place to do this.
    }
    
    func testBC() {
        // boundary conditions should be that psi(0) = 0 = psi(L)
        var solver = SchrodingerSolver()
        let xMin = 0.0, xMax = 1.0, amplitude = 0.0
        let potential : PotentialType = .squareWell, potentialList = getPotential(xMin: xMin, xMax: xMax, steps: 250, choice: potential, amplitude: amplitude)
        solver.boundaryValProblem(a: xMax - xMin, steps: 250, Vf: potentialList, ic: (psi: 0, psip: 1), eMin: 0.5, eMax: 5, eStride: 0.75)
        XCTAssertEqual(solver.totalFuncToPlot[0][0][1], 0.0, accuracy: 1e-10)
        XCTAssertEqual(solver.totalFuncToPlot[0].last[1], 0.0, accuracy: 1e-10)
    }

    override func tearDownWithError() throws {
        // Put teardown code here. This method is called after the invocation of each test method in the class.
    }

    func testExample() throws {
        // UI tests must launch the application that they test.
        let app = XCUIApplication()
        app.launch()

        // Use recording to get started writing UI tests.
        // Use XCTAssert and related functions to verify your tests produce the correct results.
    }

    func testLaunchPerformance() throws {
        if #available(macOS 10.15, iOS 13.0, tvOS 13.0, watchOS 7.0, *) {
            // This measures how long it takes to launch your application.
            measure(metrics: [XCTApplicationLaunchMetric()]) {
                XCUIApplication().launch()
            }
        }
    }
}
