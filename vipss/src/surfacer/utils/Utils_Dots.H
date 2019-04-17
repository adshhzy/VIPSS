#ifndef SDDOTS_DEFND
#define SDDOTS_DEFND

#include <utils/Mesh_Array.H>
#include <utils/Rn_Defs.H>
#include <utils/Utils_Color.H>

/** \class UTILSDot Utils_Dots.H utils/Utils_Dots.H
  \ingroup UtilitiesGeom
  \brief The dot, used in the dot array.

    Dot has a name (default is 1 more than last allocated, starting at
    0) a group id (used in the faces project to indicate location
    relative to lips and eyes) a color, a neutral location, and a
    current offset (used in the face project to offset the dot at each
    frame).

    Don't confuse the name with the id.

    */
class UTILSDot {

private:
   static int     m_siNames;
   int            m_iName;
   UTILSColor     m_oColor;
   int            m_oiGroupId;
   R3Pt           m_ptNeutralLoc;
   R3Vec          m_vecCurOffset;

public:
    /**@name Accessors */
    //@{
    /// 
    int                   Name( )  const { return m_iName; }
    ///
    const UTILSColor    & Color( ) const { return m_oColor; }
    ///
          UTILSColor    & Color( )       { return m_oColor; }
    ///
    const R3Pt          & Loc() const { return m_ptNeutralLoc; }
    ///
    const R3Vec         & Offset() const { return m_vecCurOffset; }

    ///
    void                 Set_name( int in_iI ) { m_iName = in_iI; }
    ///
    void                 Set_offset( const R3Vec & v ) { m_vecCurOffset = v; }
    ///
    void                 Set_point( const R3Pt & p ) { m_ptNeutralLoc = p; }
    ///
    void                 Set_point( const R2Pt & p ) { m_ptNeutralLoc = R3Pt( p[0], p[1], 0 ); }
    ///
    void                 Set_groupid( int id ) { m_oiGroupId = id; }
    ///
    int                  GroupId( ) const { return m_oiGroupId; }
    //@}

    /**@name Comparison operators for group id 

       These routines are used in the facial code and examine the
       group id to determine the type of dot. 

       Dots are compatable if they're not one above and one below a
       feature, where the features are the eyes and lips.

    */
    //@{
    ///
    WINbool   Is_lip( ) const;
    ///
    WINbool   Is_eye( ) const;
    ///
    WINbool   Is_above_lip( ) const;
    ///
    WINbool   Is_below_lip( ) const;
    /// Pass in the group id
    WINbool   Is_compatable( int id ) const;
    /// Pas in the dot itself
    WINbool   Is_compatable( const UTILSDot & d ) const { return Is_compatable( d.m_oiGroupId ); }
    //@}
   
    /**@name Constructors and assignments */
    //@{
    ///
    UTILSDot( const UTILSDot & );
    /// Zero everything, name is default name plus one, group id is -1
    UTILSDot( ) ;
    ///
    UTILSDot & operator=( const UTILSDot & d );
    ///
   ~UTILSDot( ) {}
    //@}

    /**@name Read and write */
    //@{
    ///
    void                 Read( ifstream & in );
    /// Only read in point
    void                 Read_locs_only( ifstream & in );
    ///
    void                 Write( ofstream & ) const;
    //@}

};

/** \class UTILSDotArray Utils_Dots.H utils/Utils_Dots.H 
  \ingroup UtilitiesGeom
  \brief Drawable dots with an offset, name, and a color.

    Each dot has a name, a group id, a location, a color, and an
    offset. There are routines for drawing these with the appropriate
    names, colors, and with or without the offset. Originally designed
    to model the dots which we glued on the face, they have since
    evolved into being the drawable equivalant to point arrays.

    Drawing routines are in the utils/OpenGL library. Options
    include drawing the names, drawing the dots at their neutral
    location plus an offset, or drawing the offset as a line from the
    dot.  

      See include/utils/Ogl_DrawDots.H for drawing options.

<br>   Files: 
- include/utils/Utils_Dots.H
- include/utils/Utils_Dots_i.H
- include/utils/Ogl_DrawDots.H
- src/utils/Utils_Dots.cpp
- src/utils/OpenGL/OGLDrawDots.cpp

Note: To use the DrawDots stuff you have to link with the utils/OpenGL library.

*/
class UTILSDotArray {
public:
    /**@name Accessors */
    //@{
    ///
    int                N( ) const { return m_aoDots.num( ); }
    ///
          UTILSDot &   operator[  ]( int in_iI )       { return m_aoDots[ in_iI ]; }
    ///
    const UTILSDot &   operator[  ]( int in_iI ) const { return m_aoDots[ in_iI ]; }
    ///
    UTILSDot &         Last( ) { return m_aoDots.last( ); }
    ///
    const UTILSDot &   Last( ) const { return m_aoDots.last( ); }
    /// Turn into a list of points (does not add offset in)
    Array <R3Pt>       PointsOnly( ) const;
    ///
    pair< R3Pt, double > CenterScale() const;
    //@}

    /**@name Adding, deleting, and moving dots */
    //@{
    /// Set to as many points as needed. The dots are undefined after this
    void             Need( int iN );
    /// Add with point, color, name
    UTILSDot &       Add( const R3Pt & p, const UTILSColor &, int in_iName = -1 );
    /// Add with point and color
    UTILSDot &       Add( const R2Pt & p, const UTILSColor &, int in_iName = -1 );
    
    /// Add with point, vector, and color
    UTILSDot &       Add( const R2Pt & p, const R2Vec &v, const UTILSColor & );
    /// Add with point, vector, and color
    UTILSDot &       Add( const R3Pt & p, const R3Vec &v, const UTILSColor & );
    /// Add dot
    void             Add( const UTILSDot & d ) { m_aoDots += d; }
    /// Add dots
    void             Add( const UTILSDotArray & in_rcoDots );
    
    /// Add points, vector between them, and color, no loop
    void             Add( const Array<R2Pt> & in_apt, const UTILSColor & );
    /// Add points, vector between them, and color, no loop
    void             Add( const Array<R2Pt> & in_apt, const Array < UTILSColor > & in_arcoColor );
    /// Add points, vector between them, and color, no loop
    void             Add( const Array<R2Pt> &in_apt, const Array < double > & in_adColVals,
                              const UTILSColor::RampType in_type = UTILSColor::INCREASING_RGB );
    /// Add points, vector between them, and color, no loop
    void             Add( const Array<R2Pt> & in_apt, const UTILSColor::RampType in_type );
    /// Add points, vector between them, and color, no loop
    void             Add( const Array<R2Pt> & in_apt, const UTILSColor::Color in_color ){
                        Add( in_apt, UTILSColor(in_color) );
    };
    
    /// Add points, vector between them, and color, no loop
    void             Add( const Array<R3Pt> & in_apt, const UTILSColor & );
    /// Add points, vector between them, and color, no loop
    void             Add( const Array<R3Pt> & in_apt, const Array < UTILSColor > & in_arcoColor );
    /// Add points, vector between them, and color, no loop
    void             Add( const Array<R3Pt> & in_apt, const Array < double > & in_adColVals,
                              const UTILSColor::RampType in_type = UTILSColor::INCREASING_RGB );
    /// Add points, vector between them, and color, no loop
    void             Add( const Array<R3Pt> & in_apt, const UTILSColor::RampType in_type );
    /// Add points, vector between them, and color, no loop
    void             Add( const Array<R3Pt> & in_apt, const UTILSColor::Color in_color ){
                        Add( in_apt, UTILSColor(in_color) );
    };
    
    void             Add( const Array<R2Pt> &in_apt, const Array<R2Vec> &in_avec, const UTILSColor &in_rcoColor );
    void             Add( const Array<R2Pt> &in_apt, const Array<R2Vec> &in_avec,
                          const UTILSColor::RampType in_type = UTILSColor::INCREASING_RGB );
    void             Add( const Array<R2Pt> &in_apt, const Array<R2Vec> &in_avec, const Array< UTILSColor > & in_arcoColor );
    void             Add( const Array<R2Pt> &in_apt, const Array<R2Vec> &in_avec,
                         const Array < double > & in_adColVals, const UTILSColor::RampType in_type = UTILSColor::INCREASING_RGB );
    void             Add( const Array<R2Pt> &in_apt, const Array<R2Vec> &in_avec, const UTILSColor::Color in_color ){
                            Add( in_apt, in_avec, UTILSColor(in_color) );
    };

    void             Add( const Array<R3Pt> &in_apt, const Array<R3Vec> &in_avec, const UTILSColor &in_rcoColor );
    void             Add( const Array<R3Pt> &in_apt, const Array<R3Vec> &in_avec,
                          const UTILSColor::RampType in_type = UTILSColor::INCREASING_RGB );
    void             Add( const Array<R3Pt> &in_apt, const Array<R3Vec> &in_avec, const Array< UTILSColor > & in_arcoColor );
    void             Add( const Array<R3Pt> &in_apt, const Array<R3Vec> &in_avec,
                          const Array < double > & in_adColVals, const UTILSColor::RampType in_type = UTILSColor::INCREASING_RGB );
    void             Add( const Array<R3Pt> &in_apt, const Array<R3Vec> &in_avec, const UTILSColor::Color in_color ){
                            Add( in_apt, in_avec, UTILSColor(in_color) );
    };
    
    /// Add two points, vector between them, and color
    void             AddPair( const R2Pt & p1, const R2Pt & p2, const UTILSColor & );
    /// Add two points, vector between them, and color
    void             AddPair( const Array<R2Pt> & in_apt1, const Array<R2Pt> & in_apt2, const UTILSColor & );
    /// Add two points, vector between them, and color
    void             AddPair( const Array<R2Pt> & in_apt1, const Array<R2Pt> & in_apt2, const Array < UTILSColor > & in_arcoColor );
    /// Add two points, vector between them, and color
    void             AddPair( const Array<R2Pt> & in_apt1, const Array<R2Pt> & in_apt2, const Array < double > & in_adColVals,
                             const UTILSColor::RampType in_type = UTILSColor::INCREASING_RGB );
    /// Add two points, vector between them, and color
    void             AddPair( const Array<R2Pt> & in_apt1, const Array<R2Pt> & in_apt2, const UTILSColor::RampType in_type );
    /// Add two points, vector between them, and color
    void             AddPair( const Array<R2Pt> & in_apt1, const Array<R2Pt> & in_apt2, const UTILSColor::Color in_color ){
                            AddPair( in_apt1, in_apt2, UTILSColor(in_color) );
    };

    
    /// Add two points, vector between them, and color
    void             AddPair( const R3Pt & p1, const R3Pt & p2, const UTILSColor & );
    /// Add two points, vector between them, and color
    void             AddPair( const Array<R3Pt> & in_apt1, const Array<R3Pt> & in_apt2, const UTILSColor & );
    /// Add two points, vector between them, and color
    void             AddPair( const Array<R3Pt> & in_apt1, const Array<R3Pt> & in_apt2, const Array < UTILSColor > & in_arcoColor );
    /// Add two points, vector between them, and color
    void             AddPair( const Array<R3Pt> & in_apt1, const Array<R3Pt> & in_apt2, const Array < double > & in_adColVals,
                             const UTILSColor::RampType in_type = UTILSColor::INCREASING_RGB );
    /// Add two points, vector between them, and color
    void             AddPair( const Array<R3Pt> & in_apt1, const Array<R3Pt> & in_apt2, const UTILSColor::RampType in_type );
    /// Add two points, vector between them, and color
    void             AddPair( const Array<R3Pt> & in_apt1, const Array<R3Pt> & in_apt2, const UTILSColor::Color in_color ){
                            AddPair( in_apt1, in_apt2, UTILSColor(in_color) );
    };

    
    /// Add connected points, vector between them, and color, no loop
    void             AddCurve( const Array<R2Pt> & in_apt, const UTILSColor & );
    /// Add connected points, vector between them, and color, no loop
    void             AddCurve( const Array<R2Pt> & in_apt, const Array < UTILSColor > & in_arcoColor );
    /// Add connected points, vector between them, and color, no loop
    void             AddCurve( const Array<R2Pt> &in_apt, const Array < double > & in_adColVals,
                               const UTILSColor::RampType in_type = UTILSColor::INCREASING_RGB );
    /// Add connected points, vector between them, and color, no loop
    void             AddCurve( const Array<R2Pt> & in_apt, const UTILSColor::RampType in_type );
    /// Add connected points, vector between them, and color, no loop
    void             AddCurve( const Array<R2Pt> & in_apt, const UTILSColor::Color in_color ){
                        AddCurve( in_apt, UTILSColor(in_color) );
    };
    
    /// Add connected points, vector between them, and color, no loop
    void             AddCurve( const Array<R3Pt> & in_apt, const UTILSColor & );
    /// Add connected points, vector between them, and color, no loop
    void             AddCurve( const Array<R3Pt> & in_apt, const Array < UTILSColor > & in_arcoColor );
    /// Add connected points, vector between them, and color, no loop
    void             AddCurve( const Array<R3Pt> & in_apt, const Array < double > & in_adColVals,
                               const UTILSColor::RampType in_type = UTILSColor::INCREASING_RGB );
    /// Add connected points, vector between them, and color, no loop
    void             AddCurve( const Array<R3Pt> & in_apt, const UTILSColor::RampType in_type );
    /// Add connected points, vector between them, and color, no loop
    void             AddCurve( const Array<R3Pt> & in_apt, const UTILSColor::Color in_color ){
                        AddCurve( in_apt, UTILSColor(in_color) );
    };
    
    /// Add connected points, vector between them, and color, loop
    void             AddPolygon( const Array<R2Pt> & in_apt, const UTILSColor & );
    /// Add connected points, vector between them, and color, loop
    void             AddPolygon( const Array<R2Pt> &in_apt, const Array < UTILSColor > & in_arcoColor );
    /// Add connected points, vector between them, and color, loop
    void             AddPolygon( const Array<R2Pt> &in_apt, const Array < double > & in_adColVals,
                                 const UTILSColor::RampType in_type = UTILSColor::INCREASING_RGB );
    /// Add connected points, vector between them, and color, loop
    void             AddPolygon( const Array<R2Pt> &in_apt, const UTILSColor::RampType in_type );
    /// Add connected points, vector between them, and color, loop
    void             AddPolygon( const Array<R2Pt> & in_apt, const UTILSColor::Color in_color ){
                        AddPolygon( in_apt, UTILSColor(in_color) );
    };
    
    /// Add connected points, vector between them, and color, loop
    void             AddPolygon( const Array<R3Pt> & in_apt, const UTILSColor & );
    /// Add connected points, vector between them, and color, loop
    void             AddPolygon( const Array<R3Pt> &in_apt, const Array < UTILSColor > & in_arcoColor );
    /// Add connected points, vector between them, and color, loop
    void             AddPolygon( const Array<R3Pt> &in_apt, const Array < double > & in_adColVals,
                                 const UTILSColor::RampType in_type = UTILSColor::INCREASING_RGB );
    /// Add connected points, vector between them, and color, loop
    void             AddPolygon( const Array<R3Pt> &in_apt, const UTILSColor::RampType in_type );
    /// Add connected points, vector between them, and color, loop
    void             AddPolygon( const Array<R3Pt> & in_apt, const UTILSColor::Color in_color ){
                        AddPolygon( in_apt, UTILSColor(in_color) );
    };
    
    ///
    void             RemoveAll( ) { m_aoDots.clearcompletely( ); }
    ///
    void             RemoveAt( int which );
    /// Move all the dots
    void             Move( const R3Matrix & f );
    /// Center at origin, +- 0.5
    void             Center();
    /** Rearrange the dots so dots with same group ids are together.
        Also renumbers */
    void             Order_by_group_ids( );
    /// Renumber dots, setting name
    void             Set_names( );
    //@}

    /**@name Operations */
    //@{
    /// The average distance between dots
    double           Avg_dist( ) const;
    /// Find the closest dot
    int              Closest( const R3Pt & p ) const;
	/// Find the closest dot in screen space (located in OpenGL)
	int              ClosestScreenSpace( const R2Pt &, double in_dNear ) const;
    //@}


    /**@name Constructors and assignments */
    //@{
    ///
    UTILSDotArray( );
    ///
    UTILSDotArray( const UTILSDotArray & ds );
    ///
    UTILSDotArray & operator=( const UTILSDotArray & dots );
    ///
    ~UTILSDotArray( ) {}
    //@}

    /**@name Read, write */
    //@{
    ///
    int  Read( ifstream & in );
    ///
    void Write( ofstream & out ) const;
    /// Just the points
    void Write_locs_only( ofstream & out ) const;
    /// Opens up file, reads, and closes file
    void Read( const char * );
    /// Opens up file, writes, and closes file
    void Write( const char * ) const;
    //@}

private:

   Array <UTILSDot>        m_aoDots;

   int      Read_locs_only( ifstream & in );



};

#include "Utils_Dots_i.H"

#endif
