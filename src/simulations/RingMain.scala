package simulations

import scalafx.application.JFXApp
import scalafx.scene.canvas.Canvas
import scalafx.scene.Scene
import scalafx.animation.AnimationTimer

object RingMain extends JFXApp {
  val cellSize = 1e-4
  val massFactor = 10
  val minRadius = 2e-7
  val maxRadius = 2e-6
  def randomOrbitBody(): MutableBody = {
    val x = (math.random-0.5)*cellSize
    val y = (math.random-0.5)*cellSize
    val z = (math.random-0.5)*2*maxRadius
    val vx = (math.random-0.5)*cellSize
    val vy = -1.5*x + (math.random-0.5)*cellSize
    val radius = maxRadius*(0.1+math.random*math.random*math.random*0.9)
    val mass = massFactor * radius * radius * radius 
    new MutableBody(new MVect3(x, y, z), new MVect3(vx, vy, 0), mass, radius)
  }
  
  val sim = new NBodyMutableSim(0.00002,
      Array.fill(10000)(randomOrbitBody()))
  
  import NBodyMutableSim._
  
  val treeBuilder = (ps: Array[MutableBody]) => new GravCollTree(ps, cellSize)
  val treeAccel = (tree: GravCollTree, p: MutableBody) => tree.addAccel(p, collide)
  val bounds = (pi: MutableBody, t: Double) => slidingBrickBoundary(pi, t, cellSize, cellSize)

  stage = new JFXApp.PrimaryStage {
    title = "Particle stuff"
    scene = new Scene(800, 800) {
      val canvas = new Canvas(800, 800)
      content = canvas
      val gc = canvas.graphicsContext2D
      gc.translate(400, 400)
      gc.scale(400, -400)
      val timer: AnimationTimer = AnimationTimer { time =>
        sim.forSim(10, Seq(hillsForce), Seq(), Some(treeBuilder -> treeAccel), bounds)
        gc.clearRect(-100, -100, 200, 200)
        //        gc.strokeLine(-100,-100,100,100)
        //        gc.strokeLine(-100,100,100,-100)
        for (body <- sim.bodies) {
          gc.fillOval(2*(body.p.x - body.radius)/cellSize, 2*(body.p.y - body.radius)/cellSize, 4*body.radius/cellSize, 4*body.radius/cellSize)
        }
      }
      timer.start
    }
  }
}