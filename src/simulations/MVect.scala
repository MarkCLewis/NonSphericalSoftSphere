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
