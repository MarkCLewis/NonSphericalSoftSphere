package simulations

class NBodyMutableSim(val dt: Double, val bodies: Array[MutableBody]) {

  private var time = 0.0

  def printBodies() {
    for (b <- bodies) {
      println(s"${NBodyMutableSim.b1} $time ${b.p.x} ${b.p.y} ${b.p.z} ${b.v.x} ${b.v.y} ${b.v.z}")
    }
  }

  def forSim[T](steps: Int, singleBodyForces: Seq[MutableBody => Unit],
                pairBodyForces: Seq[(MutableBody, MutableBody) => Unit],
                treeInfo: Option[(Array[MutableBody] => T, (T, MutableBody) => Unit)] = None, // None means no tree 
                boundaryCondition: (MutableBody, Double) => Unit = (_, _) => {}): Unit = {

    for (_ <- 1 to steps) {
      treeInfo.foreach {
        case (builder, accelerator) =>
          val start = System.nanoTime()
          val tree = builder(bodies)
          println("Tree built in " + (System.nanoTime() - start) * 1e-9)
          val start2 = System.nanoTime()
          for (b <- bodies.par) accelerator(tree, b)
          println("Tree Used in " + (System.nanoTime() - start2) * 1e-9)
      }

      for (b <- bodies.par) {
        for (sbf <- singleBodyForces) sbf(b)
      }

      if (pairBodyForces.nonEmpty) {
        for {
          i <- bodies.indices
          j <- i + 1 until bodies.length
        } {
          for (pbf <- pairBodyForces) pbf(bodies(i), bodies(j))
        }
      }

      for (i <- bodies.indices) {
        bodies(i).kickStep(dt)
      }

      time += dt
      for (i <- bodies.indices) boundaryCondition(bodies(i), time)
    }
  }

}

object NBodyMutableSim {
  val cutoff = 0.8
  val sd = new MVect3(1, 0, 0) // softness direction
  val soft = 0.9
  val dentFactor = 0.8
  val b1 = 100000.0 //100.0
  val damping = 10000.0

  def bodiesByPositions(bods: Array[(Double, Double)]): Array[MutableBody] = {
    for (b <- bods) yield {
      new MutableBody(new MVect3(b._1, b._2, 0), new MVect3(0, 0, 0), 1e-10, 1)
    }
  }

  def randomBodies(num: Int): Array[MutableBody] = {
    Array.fill(num) {
      val r = math.random * num * 5 + 2
      val theta = (math.random - 0.5) * 1
      new MutableBody(new MVect3(r * math.sin(theta), r * math.cos(theta), 0), new MVect3(0, 0, 0), 1e-10, 1)
    }
  }

  def gridBodies(): Array[MutableBody] = {
    (Array.tabulate(7, 7) { (i, j) =>
      new MutableBody(new MVect3(i * 3 + math.random * 0.5, j * 3 + 2, 0), new MVect3(0, 0, 0), 1e-10, 1)
    }).flatten
  }

  def gravityForce(pi: MutableBody, pj: MutableBody) {
    val dx = pi.p.x - pj.p.x
    val dy = pi.p.y - pj.p.y
    val dz = pi.p.z - pj.p.z
    val dist = math.sqrt(dx * dx + dy * dy + dz * dz)
    val magi = pj.mass / (dist * dist * dist)
    pi.a.x -= magi * dx
    pi.a.y -= magi * dy
    pi.a.z -= magi * dz
    val magj = pi.mass / (dist * dist * dist)
    pj.a.x += magj * dx
    pj.a.y += magj * dy
    pj.a.z += magj * dz
  }

  def constantAccel(pi: MutableBody, g: MVect3) {
    pi.a.x += g.x
    pi.a.y += g.y
    pi.a.z += g.z
  }

  def planeBounce(pi: MutableBody, n: MVect3, d: Double) {
    val dist = (n dot pi.p) - d
    if (dist < pi.radius) {
      val overlap = pi.radius - dist
      val v = -(n dot pi.v)
      val mag = b1 * overlap + v * 10
      //val mag = 50 * overlap + v * 20
      pi.a.x += n.x * mag
      pi.a.y += n.y * mag
      pi.a.z += n.z * mag
    }
  }

  val n = 1
  val kappa = 1
  val n_z = 1

  def hillsForce(pi: MutableBody): Unit = {
    pi.a.x += 2 * n * pi.v.y - (kappa * kappa - 4 * n * n) * pi.p.x
    pi.a.y += -2 * pi.v.x
    pi.a.z += -n_z * n_z * pi.p.z
  }

  def slidingBrickBoundary(pi: MutableBody, time: Double, sx: Double, sy: Double): Unit = {
    val bx = sx * 0.5
    val by = sy * 0.5
    if (pi.p.y < -by) pi.p.y += sy
    else if (pi.p.y > by) pi.p.y -= sy
    // TODO fix x for sliding with time
    if (pi.p.x < -bx) {
      pi.p.x += sx
      pi.v.y -= 1.5 * sx
    } else if (pi.p.x > bx) {
      pi.p.x -= sx
      pi.v.y += 1.5 * sx
    }
  }
}