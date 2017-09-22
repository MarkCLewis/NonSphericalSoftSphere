package simulations

import scalafx.application.JFXApp
import scalafx.scene.Scene
import scalafx.scene.canvas.Canvas
import scalafx.animation.AnimationTimer

//Tests softness of different x, y orientations of the soft spot
//object Main extends App {
//
//  val vals = for (t <- 0.0 until math.Pi by 0.01) yield {
//    (math.cos(t), math.sin(t))
//  }
//
//  for ((x, y) <- vals) {
//    var sim = new NBodyMutableClass(0.0001)
//    sim.setSD(x, y)
//    var flag = true
//    while (flag) {
//      val maxAccel = sim.forSim(100, sim.softCollide)
//      //val maxAccel = sim.forSim(100,sim.collide)
//      //println(maxAccel)
//      if (maxAccel < 1e-4) {
//        flag = false
//        sim.printBodies()
//      }
//    }
//  }
//
//}

//object Main extends App {
//  val values = Array(1000, 500, 400, 300, 250, 200, 150, 100, 50, 25, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1)
//  //val values = Array(100)
//  for (v <- values) {
//    var sim = new NBodyMutableClass(0.0001, v) // v is b1
//    var flag = true
//    while (flag) {
//      //for (_ <- 1 to 1000) {
//      val maxAccel = sim.forSim(100, sim.softCollide)
//      //println(maxAccel)
//      if (maxAccel < 1e-4) {
//        flag = false
//        sim.printBodies()
//      }
//    }
//  }
//}

object Main extends JFXApp {
  val sim = new NBodyMutableSim(0.0001,
      Array(new MutableBody(new MVect3(0, 5, 0), new MVect3(0, 0, 0), 1e-10, 1),
            new MutableBody(new MVect3(0, 1, 0), new MVect3(0, 0, 0), 1e-10, 1)))
  //sim.stackBodies()
  //sim.randomBodies(10)
  
  import NBodyMutableSim._
  
  val constAccel = (pi: MutableBody) => constantAccel(pi, new MVect3(0, -1, 0))
  
          //        val p1 = new MVect3(1 / sqrt(2), 1 / sqrt(2), 0)
        //        val p2 = new MVect3(1 / (-sqrt(2)), 1 / sqrt(2), 0)
        //        planeBounce(i, p1, 0)
        //        planeBounce(i, p2, 0)

  
  val p = new MVect3(0, 1, 0)
  val plane1 = (pi: MutableBody) => planeBounce(pi, p, 0)
  val plane2 = (pi: MutableBody) => planeBounce(pi, new MVect3(1, 0, 0), -50)
  val plane3 = (pi: MutableBody) => planeBounce(pi, new MVect3(-1, 0, 0), -50)


  stage = new JFXApp.PrimaryStage {
    title = "Particle stuff"
    scene = new Scene(800, 800) {
      val canvas = new Canvas(800, 800)
      content = canvas
      val gc = canvas.graphicsContext2D
      gc.translate(400, 400)
      gc.scale(5, -5)
      val timer: AnimationTimer = AnimationTimer { time =>
        sim.forSim(100, Seq(constAccel, plane1, plane2, plane3), Seq(warpedCollide))
        gc.clearRect(-100, -100, 200, 200)
        //        gc.strokeLine(-100,-100,100,100)
        //        gc.strokeLine(-100,100,100,-100)
        for (body <- sim.bodies) {
          gc.fillOval(body.p.x - 1, body.p.y - 1, 2, 2)
        }
      }
      timer.start
    }
  }
}


// sd = softness direction, cutoff, and soft = softness factor

// this if the value that will replace b1... "soft b1"
// if (sd dot n > cutoff) {
// b1 * ( 1 - soft * f(sd))
// } else { b1 }

// f(x) = (x - cutoff) / (1 - cutoff)

/* 
 * Timestep = 0.000001
 * steps = 10000
 * 
 * b1 = 1000
 * 6.290000000525299 0.0 2.9980008704216035 0.0 0.0 -1.2451621835194464E-5 0.0
 * 6.290000000525299 0.0 0.9979988687751139 0.0 0.0 -2.9970817137476443E-5 0.0
 * 
 * b1 = 500
 * 6.260000000521106 0.0 2.996003570237496 0.0 0.0 -3.135273724638581E-4 0.0
 * 6.260000000521106 0.0 0.9960002310429089 0.0 0.0 2.792753758878369E-4 0.0
 * 
 * b1 = 100
 * 4.920000000333803 0.0 2.9800002929436498 0.0 0.0 9.353180473816068E-5 0.0
 * 4.920000000333803 0.0 0.9800002137323449 0.0 0.0 -9.201680930065943E-5 0.0
 * 
 * b1 = 50
 * -- not possible -- gains too much energy at 50 ???
 */