package simulations
import Math.sqrt

class MVect3(var x: Double, var y: Double, var z: Double) {

  def zero(): Unit = {
    x = 0.0
    y = 0.0
    z = 0.0
  }

  def dot(m: MVect3): Double = {
    x * m.x + y * m.y + z * m.z
  }

  def normalize() {
    val mag = sqrt(this dot this)
    x /= mag
    y /= mag
    z /= mag
  }

  def -(m: MVect3): MVect3 = new MVect3(x - m.x, y - m.y, z - m.z)

  def mag() = {
    sqrt(this dot this)
  }
}

class MutableBody(
  val p: MVect3,
  val v: MVect3,
  val mass: Double,
  val radius: Double)

class NBodyMutableClass(val dt: Double, val bodies: Array[MutableBody]) {
  var cutoff = 0.8
  var sd = new MVect3(1, 0, 0) // softness direction
  var soft = 0.9
  var b1 = 100.0;

  private val accel = Array.fill(bodies.length)(new MVect3(0, 0, 0))
  private var time = 0.0

  def setB(b: Double) {
    b1 = b
  }

  def bodiesByPositions(bods: Array[(Double, Double)]): Array[MutableBody] = {
    for (b <- bods) yield {
      new MutableBody(new MVect3(b._1, b._2, 0), new MVect3(0, 0, 0), 1e-10, 1)
    }
  }

  def setSD(x: Double, y: Double) {
    sd = new MVect3(x, y, 0)
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

  //  initBodies()
  //  def initBodies(): Unit = {
  //    bodies(0).mass = 1.0
  //    for (i <- 1 until numBodies) {
  //      bodies(i).p.x = i
  //      bodies(i).v.y = math.sqrt(1.0 / i)
  //    }
  //  }

  def printBodies() {
    for (b <- bodies) {
      println(s"$b1 $time ${b.p.x} ${b.p.y} ${b.p.z} ${b.v.x} ${b.v.y} ${b.v.z}")
    }
  }

  def forSim(steps: Int, singleBodyForces: Seq[Int => Unit], pairBodyForces: Seq[(Int, Int) => Unit]): Double = {
    var maxAccel = 0.0
    for (_ <- 1 to steps) {

      for (i <- bodies.indices) {
        accel(i).zero()
        for (sbf <- singleBodyForces) sbf(i)
      }

      for {
        i <- bodies.indices
        j <- i + 1 until bodies.length
      } {
        for (pbf <- pairBodyForces) pbf(i, j)
      }
      for (i <- bodies.indices) {
        val p = bodies(i)
        val a = accel(i).mag()
        if (a > maxAccel) maxAccel = a
        p.v.x += accel(i).x * dt
        p.v.y += accel(i).y * dt
        p.v.z += accel(i).z * dt
        p.p.x += p.v.x * dt
        p.p.y += p.v.y * dt
        p.p.z += p.v.z * dt
      }
      time += dt
    }
    maxAccel
  }

  def f(x: Double): Double = (x - cutoff) / (1 - cutoff) // fixxxxxxxxxxx

  def softCollide(i: Int, j: Int) {
    val pi = bodies(i)
    val pj = bodies(j)
    val dx = pi.p.x - pj.p.x
    val dy = pi.p.y - pj.p.y
    val dz = pi.p.z - pj.p.z
    val dist = math.sqrt(dx * dx + dy * dy + dz * dz)
    if (dist < (pi.radius + pj.radius)) {
      val overlap = (pi.radius + pj.radius) - dist
      val n = new MVect3(dx, dy, dz)
      n.normalize
      val softFactor = sd dot n
      val softB = if (softFactor > cutoff) b1 * (1 - soft * f(softFactor)) else b1

      val v = -(n dot (pi.v - pj.v))
      val mag = -softB * overlap - v * 10
      accel(i).x -= n.x * mag //need to divide by mass
      accel(i).y -= n.y * mag
      accel(i).z -= n.z * mag
      accel(j).x += n.x * mag
      accel(j).y += n.y * mag
      accel(j).z += n.z * mag
    }
  }

  def collide(i: Int, j: Int) {
    val pi = bodies(i)
    val pj = bodies(j)
    val dx = pi.p.x - pj.p.x
    val dy = pi.p.y - pj.p.y
    val dz = pi.p.z - pj.p.z
    val dist = math.sqrt(dx * dx + dy * dy + dz * dz)
    if (dist < (pi.radius + pj.radius)) {
      val overlap = (pi.radius + pj.radius) - dist
      val n = new MVect3(dx, dy, dz)
      n.normalize
      val v = -(n dot (pi.v - pj.v))
      val mag = -b1 * overlap - v * 10
      accel(i).x -= n.x * mag //need to divide by mass
      accel(i).y -= n.y * mag
      accel(i).z -= n.z * mag
      accel(j).x += n.x * mag
      accel(j).y += n.y * mag
      accel(j).z += n.z * mag
    }
  }

  def gravityForce(i: Int, j: Int) {
    val pi = bodies(i)
    val pj = bodies(j)
    val dx = pi.p.x - pj.p.x
    val dy = pi.p.y - pj.p.y
    val dz = pi.p.z - pj.p.z
    val dist = math.sqrt(dx * dx + dy * dy + dz * dz)
    val magi = pj.mass / (dist * dist * dist)
    accel(i).x -= magi * dx
    accel(i).y -= magi * dy
    accel(i).z -= magi * dz
    val magj = pi.mass / (dist * dist * dist)
    accel(j).x += magj * dx
    accel(j).y += magj * dy
    accel(j).z += magj * dz
  }

  def constantAccel(i: Int, g: MVect3) {
    accel(i).x += g.x
    accel(i).y += g.y
    accel(i).z += g.z
  }

  def planeBounce(i: Int, n: MVect3, d: Double) {
    val pi = bodies(i)
    val dist = (n dot pi.p) - d
    if (dist < pi.radius) {
      val overlap = pi.radius - dist
      val v = -(n dot pi.v)
      val mag = b1 * overlap + v * 10
      //val mag = 50 * overlap + v * 20
      accel(i).x += n.x * mag
      accel(i).y += n.y * mag
      accel(i).z += n.z * mag
    }
  }

  val n = 1
  val kappa = 1
  val n_z = 1

  def hillsForce(i: Int): Unit = {
    val pi = bodies(i)
    accel(i).x += 2 * n * pi.v.y - (kappa * kappa - 4 * n * n) * pi.p.x
    accel(i).y += -2 * pi.v.x
    accel(i).z += -n_z * n_z * pi.p.z
  }
}