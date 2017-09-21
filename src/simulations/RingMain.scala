package simulations

import scalafx.application.JFXApp
import scalafx.scene.canvas.Canvas
import scalafx.scene.Scene
import scalafx.animation.AnimationTimer

object RingMain extends JFXApp {
  val cellWidth = 1e-4
  val cellHeight = 1e-4
  val pmass = 1e-17
  val pradius = 1e-6
  def randomOrbitBody(): MutableBody = {
    val x = (math.random-0.5)*cellWidth
    val y = (math.random-0.5)*cellHeight
    val z = (math.random-0.5)*2*pradius
    val vx = (math.random-0.5)*cellWidth
    val vy = -1.5*x + (math.random-0.5)*cellWidth
    new MutableBody(new MVect3(x, y, z), new MVect3(vx, vy, 0), pmass, pradius)
  }
  
  val sim = new NBodyMutableClass(0.00001,
      Array.fill(1000)(randomOrbitBody()))
  
  val bounds = (i: Int, t: Double) => sim.slidingBrickBoundary(i, t, cellWidth, cellHeight)

  stage = new JFXApp.PrimaryStage {
    title = "Particle stuff"
    scene = new Scene(800, 800) {
      val canvas = new Canvas(800, 800)
      content = canvas
      val gc = canvas.graphicsContext2D
      gc.translate(400, 400)
      gc.scale(400, -400)
      val timer: AnimationTimer = AnimationTimer { time =>
        sim.forSim(10, Seq(sim.hillsForce), Seq(sim.collide, sim.gravityForce), bounds)
        gc.clearRect(-100, -100, 200, 200)
        //        gc.strokeLine(-100,-100,100,100)
        //        gc.strokeLine(-100,100,100,-100)
        for (body <- sim.bodies) {
          gc.fillOval(2*(body.p.x - body.radius)/cellWidth, 2*(body.p.y - body.radius)/cellHeight, 4*body.radius/cellWidth, 4*body.radius/cellHeight)
        }
      }
      timer.start
    }
  }
}