package simulations

class MutableBody(
    val p: MVect3,
    val v: MVect3,
    val mass: Double,
    val radius: Double,
    val a: MVect3 = new MVect3(0, 0, 0)) {

  def kickStep(dt: Double): Unit = {
        v.x += a.x * dt
        v.y += a.y * dt
        v.z += a.z * dt
        p.x += v.x * dt
        p.y += v.y * dt
        p.z += v.z * dt
        a.x = 0.0
        a.y = 0.0
        a.z = 0.0
  }
}
