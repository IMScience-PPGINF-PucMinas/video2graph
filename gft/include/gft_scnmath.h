
#ifndef _GFT_SCNMATH_H_
#define _GFT_SCNMATH_H_

#include "gft_common.h"
#include "gft_scene.h"

namespace gft{
  namespace Scene32{

    sScene32 *Sub(sScene32 *scn1, sScene32 *scn2);
    /**
     * Inplace version.
     */
    void    Subinplace(sScene32 *scn1, sScene32 *scn2);

    sScene32 *Add(sScene32 *scn1, sScene32 *scn2);
    sScene32 *Add(sScene32 *scn, int value);
    /**
     * Inplace version.
     */
    void     Addinplace(sScene32 *scn1, sScene32 *scn2);

    sScene32 *Mult(sScene32 *scn1, sScene32 *scn2);

    sScene32 *Or(sScene32 *scn1, sScene32 *scn2);
    /**
     * Inplace version.
     */
    void     Orinplace(sScene32 *scn1, sScene32 *scn2);
    sScene32 *And(sScene32 *scn1, sScene32 *scn2);
    sScene32 *XOr(sScene32 *scn1, sScene32 *scn2);

    sScene32 *Complement(sScene32 *scn);
    sScene32 *Abs(sScene32 *scn);
    /**
     * Inplace version.
     */
    void     Negateinplace(sScene32 *scn);

  } //end Scene32 namespace


  namespace Scene16{

    sScene16 *Sub(sScene16 *scn1, sScene16 *scn2);
    /**
     * Inplace version.
     */
    void     Subinplace(sScene16 *scn1, sScene16 *scn2);

    sScene16 *Add(sScene16 *scn1, sScene16 *scn2);
    sScene16 *Add(sScene16 *scn, ushort value);
    /**
     * Inplace version.
     */
    void     Addinplace(sScene16 *scn1, sScene16 *scn2);

    sScene16 *Mult(sScene16 *scn1, sScene16 *scn2);

    sScene16 *Or(sScene16 *scn1, sScene16 *scn2);
    /**
     * Inplace version.
     */
    void     Orinplace(sScene16 *scn1, sScene16 *scn2);
    sScene16 *And(sScene16 *scn1, sScene16 *scn2);
    sScene16 *XOr(sScene16 *scn1, sScene16 *scn2);

    sScene16 *Complement(sScene16 *scn);

  } //end Scene16 namespace


  namespace Scene8{

    sScene8 *Sub(sScene8 *scn1, sScene8 *scn2);
    /**
     * Inplace version.
     */
    void    Subinplace(sScene8 *scn1, sScene8 *scn2);

    sScene8 *Add(sScene8 *scn1, sScene8 *scn2);
    sScene8 *Add(sScene8 *scn, uchar value);
    /**
     * Inplace version.
     */
    void    Addinplace(sScene8 *scn1, sScene8 *scn2);

    sScene8 *Mult(sScene8 *scn1, sScene8 *scn2);

    sScene8 *Or(sScene8 *scn1, sScene8 *scn2);
    /**
     * Inplace version.
     */
    void    Orinplace(sScene8 *scn1, sScene8 *scn2);
    sScene8 *And(sScene8 *scn1, sScene8 *scn2);
    sScene8 *XOr(sScene8 *scn1, sScene8 *scn2);

    sScene8 *Complement(sScene8 *scn);

  } //end Scene8 namespace


  namespace Scene{

    sScene *Sub(sScene *scn1, sScene *scn2);
    /**
     * Inplace version.
     */
    void   Subinplace(sScene *scn1, sScene *scn2);

    sScene *Add(sScene *scn1, sScene *scn2);
    sScene *Add(sScene *scn, int value);
    /**
     * Inplace version.
     */
    void   Addinplace(sScene *scn1, sScene *scn2);

    sScene *Mult(sScene *scn1, sScene *scn2);

    sScene *Or(sScene *scn1, sScene *scn2);
    /**
     * Inplace version.
     */
    void   Orinplace(sScene *scn1, sScene *scn2);
    sScene *And(sScene *scn1, sScene *scn2);
    sScene *XOr(sScene *scn1, sScene *scn2);

    sScene *Complement(sScene *scn);

  } //end Scene namespace


} //end gft namespace

#endif

