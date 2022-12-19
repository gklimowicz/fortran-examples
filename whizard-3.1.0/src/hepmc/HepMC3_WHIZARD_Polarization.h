#include <iostream>
#include <cmath>

static const double HepMC_pi = 3.14159265358979323846;  // copy of pi from CLHEP

class Polarization {
  
  /// print polarization information
  friend std::ostream& operator<<( std::ostream&, const Polarization& );
  
 public:
  /// default constructor
  Polarization( );
  /// constructor requiring at least one value
  Polarization( double theta, double phi = 0 );
  /// construct from another polarization object
  Polarization( const Polarization& inpolar );
  virtual       ~Polarization() {}
  
  /// swap
  void swap( Polarization & other);
  /// make a copy
  Polarization& operator=( const Polarization& inpolar );
  /// equality requires that theta and phi are equal
  bool          operator==( const Polarization& ) const;
  /// inequality results if either theta or phi differ
  bool          operator!=( const Polarization& ) const;
  
  /// print theta and phi
  void          print( std::ostream& ostr = std::cout ) const;
  
  ////////////////////
  // access methods //
  ////////////////////
  double        theta() const;    //!< returns polar angle in radians
  double        phi() const;      //!< returns azimuthal angle in radians
  bool          is_defined() const;   //!< returns true if the Polarization has been defined
  
  /// set polar angle in radians 
  double        set_theta( double theta );
  /// set azimuthal angle in radians 
  double        set_phi( double phi );
  /// set both polar and azimuthal angles in radians
  void          set_theta_phi( double theta, double phi );
  /// declares the Polarization as undefined and zeros the values
  void          set_undefined();
  
 private:
  /// private method to return a polar angle in the correct range
  double valid_theta( double theta );
  /// private method to return an azimuthal angle in the correct range
  double valid_phi( double phi );
  
 private:
  double m_theta; //polar angle of polarization in radians 0< theta <pi
  double m_phi;   //azimuthal angle of polarization in rad. 0< phi <2pi
  bool   m_defined; //used to flag if the Polarization has been defined
};

///////////////////////////
// INLINE Access Methods //
///////////////////////////

inline double Polarization::theta() const { return m_theta; }
inline double Polarization::phi() const { return m_phi; }

///////////////////////////
// INLINE Operators      //
///////////////////////////

inline bool Polarization::operator==( const Polarization& a ) const 
{
  return ( a.theta() == this->theta() && a.phi() == this->phi() && a.is_defined() == this->is_defined() );
}

inline bool Polarization::operator!=(const Polarization& a ) const 
{
  return !( a == *this );
}
