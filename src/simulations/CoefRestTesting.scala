package simulations

/**
 * The purpose of this code is to test the coefficient of restitution and related values for simple systems where
 * two bodies collide with one another. We want to explore a high dimensional parameter space that includes:
 *  - b1: constant determining particle hardness
 *  - damping: constant determining the damping of particle velocities during overlap
 *  - dt: time step
 *  - impact velocity
 *
 * For the non-spherical collisions we also have the cutoff and the softness/dent factors.
 */
object CoefRestTesting extends App {

  softSpheres()
  
  def plainSpheres(): Unit = {
    val b1s = Array(100.0, 300.0, 1000.0, 3000.0, 10000.0, 30000.0, 100000.0, 300000.0, 1000000.0)
    val dampings = b1s.map(_ / 10)
    val vels = Array(1.0, 3.0, 10.0, 30.0, 100.0)
    val dts = Array(1e-3, 3e-4, 1e-4, 3e-5, 1e-5)

    val pw = new java.io.PrintWriter("coefRestSphere.txt")
    for (b1 <- b1s; damping <- dampings; vel <- vels; dt <- dts) {
      val collForce = new SphereCollide(b1, damping)
      val (ceofRest, pdepth, stick, passThru) = runTest(collForce, vel, dt, 1.0)
      pw.println(s"$b1 $damping $vel $dt $ceofRest $pdepth $stick $passThru")
      println(s"$b1\t$damping\t$vel\t$dt\t$ceofRest\t$pdepth\t$stick\t$passThru")
    }
    pw.close
  }

  def softSpheres(): Unit = {
    val b1s = Array(1000.0, 3000.0, 10000.0, 30000.0, 100000.0)
    val dampings = b1s.map(_ / 10)
    val vels = Array(1.0, 3.0, 10.0, 30.0)
    val dts = Array(1e-4, 3e-5, 1e-5)
    val softenings = Array(0.01, 0.03, 0.1, 0.3, 0.5, 0.8)

    val pw = new java.io.PrintWriter("coefRestSoft.txt")
    for (b1 <- b1s; damping <- dampings; vel <- vels; dt <- dts; soft <- softenings) {
      val collForce = new SoftenedSphereCollide(b1, damping, Array(SofteningData(new MVect3(1, 0, 0), 0.8, soft),
          SofteningData(new MVect3(-1, 0, 0), 0.8, soft)))
      val (ceofRest, pdepth, stick, passThru) = runTest(collForce, vel, dt, 1.0)
      pw.println(s"$b1 $damping $vel $dt $soft $ceofRest $pdepth $stick $passThru")
      println(s"$b1\t$damping\t$vel\t$dt\t$soft\t$ceofRest\t$pdepth\t$stick\t$passThru")
    }
    pw.close
  }

  def warpedSpheres(): Unit = {
    val b1s = Array(1000.0, 3000.0, 10000.0, 30000.0, 100000.0)
    val dampings = b1s.map(_ / 10)
    val vels = Array(1.0, 3.0, 10.0, 30.0)
    val dts = Array(1e-4, 3e-5, 1e-5)
    val warps = Array(0.1, 0.3, 0.5, 0.8)

    val pw = new java.io.PrintWriter("coefRestWarpedtxt")
    for (b1 <- b1s; damping <- dampings; vel <- vels; dt <- dts; warp <- warps) {
      val collForce = new WarpedCollide(b1, damping, Array(WarpingData(new MVect3(1, 0, 0), 0.8, warp),
          WarpingData(new MVect3(-1, 0, 0), 0.8, warp)))
      val (ceofRest, pdepth, stick, passThru) = runTest(collForce, vel, dt, 1.0)
      pw.println(s"$b1 $damping $vel $dt $warp $ceofRest $pdepth $stick $passThru")
      println(s"$b1\t$damping\t$vel\t$dt\t$warp\t$ceofRest\t$pdepth\t$stick\t$passThru")
    }
    pw.close
  }

  // Returns the coefficient of restitution, the maximum penetration depth, stick, pass through
  def runTest(collForce: CollisionForce, vel: Double, dt: Double, radius: Double): (Double, Double, Int, Int) = {
    val p1 = new MutableBody(new MVect3(-1.1 * radius, 0.0, 0.0), new MVect3(vel, 0.0, 0.0), 1.0, radius)
    val p2 = new MutableBody(new MVect3(1.1 * radius, 0.0, 0.0), new MVect3(-vel, 0.0, 0.0), 1.0, radius)
    val sim = new NBodyMutableSim(dt, Array(p1, p2))
    var state = 0
    var minDist = 100.0
    var stallCount = 0
    val collide = Array(collForce)
    while (state < 2) {
      sim.forSim(1, Nil, collide, None)
      val dist = (p1.p - p2.p).mag
      if (state == 0) {
        if (dist < 2 * radius) state = 1
      }
      if (state == 1) {
        minDist = minDist min dist
        if (dist > 2 * radius) state = 2
        if ((p1.v - p2.v).mag < 1e-5) {
          stallCount += 1
          if (stallCount > 100) state = 3
        }
      }
    }
    ((p1.v - p2.v).mag / (2 * vel), minDist / (2 * radius), state, if (p1.v.x > 0.0 && state != 3) 1 else 0)
  }
}