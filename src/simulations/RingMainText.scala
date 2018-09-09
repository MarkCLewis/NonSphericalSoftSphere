package simulations

import scala.concurrent.Future

object RingMainText extends App {
  val cellSize = 3e-4
  val massFactor = 10
  val minRadius = 2e-7
  val maxRadius = 2e-6
  def randomOrbitBody(): MutableBody = {
    val x = (math.random - 0.5) * cellSize
    val y = (math.random - 0.5) * cellSize
    val z = (math.random - 0.5) * 2 * maxRadius
    val vx = (math.random - 0.5) * cellSize
    val vy = -1.5 * x + (math.random - 0.5) * cellSize
    val radius = maxRadius * (0.1 + math.random * math.random * math.random * 0.9)
    val mass = massFactor * radius * radius * radius
    new MutableBody(new MVect3(x, y, z), new MVect3(vx, vy, 0), mass, radius)
  }

  val sim = new NBodyMutableSim(0.00002,
    Array.fill(100000)(randomOrbitBody()))

  import NBodyMutableSim._

  val treeBuilder = (ps: Array[MutableBody]) => new GravCollTree(ps, cellSize)
  val v1 = new MVect3(0,1,0)
  val v2 = new MVect3(0,-1,0)
  val v3 = new MVect3(1,0,0)
  val v4 = new MVect3(-1,0,0)
  v2.normalize()
  val treeAccel = (tree: GravCollTree, p: MutableBody) => {
    val collide = new OptimizedWarpedCollide(100000, 10000, Array(
	WarpingData(v1, 0.8, 0.2),
	WarpingData(v2, 0.8, 0.2),
	WarpingData(v3, 0.8, 0.2),
	WarpingData(v4, 0.8, 0.2)
	))
    tree.addAccel(p, collide, 0.0, 0.0)
    tree.addAccel(p, collide, 0.0, cellSize)
    tree.addAccel(p, collide, 0.0, -cellSize)
  }
  val bounds = (pi: MutableBody, t: Double) => slidingBrickBoundary(pi, t, cellSize, cellSize)

  for (i <- 0 to 100000) {
    sim.forSim(1000, Seq(hillsForce), Seq(), Some(treeBuilder -> treeAccel), bounds)
		sim.printBodies()
  }
}
