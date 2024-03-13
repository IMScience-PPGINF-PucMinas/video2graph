
#include "gft_adjrel.h"

namespace gft
{
  namespace AdjRel
  {

    sAdjRel *Create(int n)
    {
      sAdjRel *A = NULL;

      A = (sAdjRel *)calloc(1, sizeof(sAdjRel));
      if (A != NULL)
      {
        A->dx = gft::AllocIntArray(n);
        A->dy = gft::AllocIntArray(n);
        A->dz = gft::AllocIntArray(n);
        A->n = n;
      }
      else
      {
        gft::Error((char *)MSG1, (char *)"AdjRel::Create");
      }

      return (A);
    }

    void Destroy(sAdjRel **A)
    {
      sAdjRel *aux;

      aux = *A;
      if (aux != NULL)
      {
        if (aux->dx != NULL)
          free(aux->dx);
        if (aux->dy != NULL)
          free(aux->dy);
        if (aux->dz != NULL)
          free(aux->dy);
        free(aux);
        *A = NULL;
      }
    }

    sAdjRel *Clone(sAdjRel *A)
    {
      sAdjRel *C;
      int i;

      C = Create(A->n);
      for (i = 0; i < A->n; i++)
      {
        C->dx[i] = A->dx[i];
        C->dy[i] = A->dy[i];
        C->dz[i] = A->dz[i];
      }
      return C;
    }

    sAdjRel *Neighborhood_4()
    { /* 4-neighborhood */
      sAdjRel *A = NULL;
      A = Create(4 + 1);
      /* place central pixel at first */
      A->dx[0] = 0;
      A->dy[0] = 0;
      A->dx[1] = 1;
      A->dy[1] = 0; /* right */
      A->dx[2] = 0;
      A->dy[2] = -1; /* top */
      A->dx[3] = -1;
      A->dy[3] = 0; /* left */
      A->dx[4] = 0;
      A->dy[4] = 1; /* bottom */
      return A;
    }

    sAdjRel *Neighborhood_6()
    { /* 4-neighborhood */
      sAdjRel *A = NULL;
      A = Create(6 + 1);
      /* place central pixel at first */
      A->dx[0] = 0;
      A->dy[0] = 0;
      A->dz[0] = 0;

      A->dx[1] = 1;
      A->dy[1] = 0;
      A->dz[1] = 0; /* right */

      A->dx[2] = 0;
      A->dy[2] = -1;
      A->dz[2] = 0; /* top */

      A->dx[3] = -1;
      A->dy[3] = 0;
      A->dz[3] = 0; /* left */

      A->dx[4] = 0;
      A->dy[4] = 1;
      A->dz[4] = 0; /* bottom */

      A->dx[5] = 0;
      A->dy[5] = 0;
      A->dz[5] = -1; // Center-back

      A->dx[6] = 0;
      A->dy[6] = 0;
      A->dz[6] = 1; // Center-front

      return A;
    }

    sAdjRel *Neighborhood_8()
    { /* 8-neighborhood */
      sAdjRel *A = NULL;
      int i, dx, dy;
      A = Create(8 + 1);
      /* place central pixel at first */
      A->dx[0] = 0;
      A->dy[0] = 0;
      i = 1;
      for (dy = -1; dy <= 1; dy++)
      {
        for (dx = -1; dx <= 1; dx++)
        {
          if ((dx != 0) || (dy != 0))
          {
            A->dx[i] = dx;
            A->dy[i] = dy;
            i++;
          }
        }
      }
      return A;
    }

    sAdjRel *Neighborhood_8_counterclockwise()
    { /* 8-neighborhood */
      sAdjRel *A = NULL;
      A = Create(8 + 1);
      /* place central pixel at first */
      A->dx[0] = 0;
      A->dy[0] = 0;
      A->dx[1] = -1;
      A->dy[1] = 0;
      A->dx[2] = -1;
      A->dy[2] = 1;
      A->dx[3] = 0;
      A->dy[3] = 1;
      A->dx[4] = 1;
      A->dy[4] = 1;
      A->dx[5] = 1;
      A->dy[5] = 0;
      A->dx[6] = 1;
      A->dy[6] = -1;
      A->dx[7] = 0;
      A->dy[7] = -1;
      A->dx[8] = -1;
      A->dy[8] = -1;
      return A;
    }

    sAdjRel *Neighborhood_8_clockwise()
    { /* 8-neighborhood */
      sAdjRel *A = NULL;
      A = Create(8 + 1);
      /* place central pixel at first */
      A->dx[0] = 0;
      A->dy[0] = 0;
      A->dx[1] = -1;
      A->dy[1] = 0;
      A->dx[2] = -1;
      A->dy[2] = -1;
      A->dx[3] = 0;
      A->dy[3] = -1;
      A->dx[4] = 1;
      A->dy[4] = -1;
      A->dx[5] = 1;
      A->dy[5] = 0;
      A->dx[6] = 1;
      A->dy[6] = 1;
      A->dx[7] = 0;
      A->dy[7] = 1;
      A->dx[8] = -1;
      A->dy[8] = 1;
      return A;
    }

    sAdjRel *Circular(float r)
    {
      sAdjRel *A = NULL;
      int i, n, dx, dy, r0, r2;

      n = 0;
      r0 = (int)r;
      r2 = (int)(r * r + 0.5);
      for (dy = -r0; dy <= r0; dy++)
        for (dx = -r0; dx <= r0; dx++)
          if (((dx * dx) + (dy * dy)) <= r2)
            n++;

      A = Create(n);
      i = 1;
      for (dy = -r0; dy <= r0; dy++)
        for (dx = -r0; dx <= r0; dx++)
          if (((dx * dx) + (dy * dy)) <= r2)
          {
            if ((dx != 0) || (dy != 0))
            {
              A->dx[i] = dx;
              A->dy[i] = dy;
              i++;
            }
          }

      /* place central pixel at first */
      A->dx[0] = 0;
      A->dy[0] = 0;

      return (A);
    }

    sAdjRel *Box(int ncols, int nrows)
    {
      sAdjRel *A = NULL;
      int i, dx, dy;

      if (ncols % 2 == 0)
        ncols++;
      if (nrows % 2 == 0)
        nrows++;

      A = Create(ncols * nrows);
      i = 1;
      for (dy = -nrows / 2; dy <= nrows / 2; dy++)
      {
        for (dx = -ncols / 2; dx <= ncols / 2; dx++)
        {
          if ((dx != 0) || (dy != 0))
          {
            A->dx[i] = dx;
            A->dy[i] = dy;
            i++;
          }
        }
      }
      /* place the central pixel at first */
      A->dx[0] = 0;
      A->dy[0] = 0;

      return (A);
    }

    int GetFrameSize(sAdjRel *A)
    {
      int sz = INT_MIN, i = 0;

      for (i = 0; i < A->n; i++)
      {
        if (abs(A->dx[i]) > sz)
          sz = abs(A->dx[i]);
        if (abs(A->dy[i]) > sz)
          sz = abs(A->dy[i]);
      }
      return (sz);
    }

    void Mult(sAdjRel *A, int val)
    {
      int i;
      for (i = 0; i < A->n; i++)
      {
        A->dx[i] *= val;
        A->dy[i] *= val;
      }
    }

    sAdjPxl *AdjPixels(sAdjRel *A, int ncols)
    {
      sAdjPxl *N;
      int i;

      N = (sAdjPxl *)calloc(1, sizeof(sAdjPxl));
      if (N != NULL)
      {
        N->dp = gft::AllocIntArray(A->n);
        N->n = A->n;
        for (i = 0; i < N->n; i++)
          N->dp[i] = A->dx[i] + ncols * A->dy[i];
      }
      else
        Error((char *)MSG1, (char *)"AdjRel::AdjPixels");
      return (N);
    }

    sAdjPxl *AdjPixels(sAdjRel *A, sImage32 *img)
    {
      return AdjPixels(A, img->ncols);
    }

    sAdjPxl *AdjPixels(sAdjRel *A, sCImage *cimg)
    {
      return AdjPixels(A, cimg->C[0]->ncols);
    }

    void DestroyAdjPxl(sAdjPxl **N)
    {
      sAdjPxl *aux;
      aux = *N;
      if (aux != NULL)
      {
        if (aux->dp != NULL)
          gft::FreeIntArray(&aux->dp);
        free(aux);
        *N = NULL;
      }
    }

    sAdjRel *RightSide(sAdjRel *A)
    {
      sAdjRel *R = NULL;
      int i;
      float d;
      /* Let p -> q be an arc represented by the increments dx,dy. Its
	 right side is given by the increments Dx = -dy/d + dx/2 and Dy =
	 dx/d + dy/2, where d=sqrt(dx²+dy²). */
      R = Create(A->n);
      for (i = 0; i < R->n; i++)
      {
        d = sqrt(A->dx[i] * A->dx[i] + A->dy[i] * A->dy[i]);
        if (d != 0)
        {
          R->dx[i] = ROUND(((float)A->dx[i] / 2.0) - ((float)A->dy[i] / d));
          R->dy[i] = ROUND(((float)A->dx[i] / d) + ((float)A->dy[i] / 2.0));
        }
      }
      return (R);
    }

    sAdjRel *LeftSide(sAdjRel *A)
    {
      sAdjRel *L = NULL;
      int i;
      float d;
      /* Let p -> q be an arc represented by the increments dx,dy. Its
	 left side is given by the increments Dx = dy/d + dx/2 and Dy =
	 -dx/d + dy/2, where d=sqrt(dx²+dy²). */
      L = Create(A->n);
      for (i = 0; i < L->n; i++)
      {
        d = sqrt(A->dx[i] * A->dx[i] + A->dy[i] * A->dy[i]);
        if (d != 0)
        {
          L->dx[i] = ROUND(((float)A->dx[i] / 2.0) + ((float)A->dy[i] / d));
          L->dy[i] = ROUND(((float)A->dy[i] / 2) - ((float)A->dx[i] / d));
        }
      }
      return (L);
    }

    sAdjRel *RightSide8(sAdjRel *A)
    {
      sAdjRel *R = NULL;
      int i;
      R = gft::AdjRel::Create(A->n);
      for (i = 0; i < R->n; i++)
      {
        if (A->dx[i] == 0 && A->dy[i] == 0)
        {
          R->dx[i] = 0;
          R->dy[i] = 0;
        }
        else if (abs(A->dx[i]) == 1 && A->dy[i] == 0)
        {
          R->dx[i] = A->dx[i];
          R->dy[i] = A->dx[i];
        }
        else if (A->dx[i] == 0 && abs(A->dy[i]) == 1)
        {
          R->dx[i] = -A->dy[i];
          R->dy[i] = A->dy[i];
        }
        else if (A->dx[i] == 1 && A->dy[i] == -1)
        {
          R->dx[i] = A->dx[i];
          R->dy[i] = 0;
        }
        else if (A->dx[i] == -1 && A->dy[i] == 1)
        {
          R->dx[i] = A->dx[i];
          R->dy[i] = 0;
        }
        else if (A->dx[i] == 1 && A->dy[i] == 1)
        {
          R->dx[i] = 0;
          R->dy[i] = A->dx[i];
        }
        else if (A->dx[i] == -1 && A->dy[i] == -1)
        {
          R->dx[i] = 0;
          R->dy[i] = A->dx[i];
        }
      }
      return (R);
    }

    sAdjRel *LeftSide8(sAdjRel *A)
    {
      sAdjRel *L = NULL;
      int i;
      L = gft::AdjRel::Create(A->n);
      for (i = 0; i < L->n; i++)
      {
        if (A->dx[i] == 0 && A->dy[i] == 0)
        {
          L->dx[i] = 0;
          L->dy[i] = 0;
        }
        else if (abs(A->dx[i]) == 1 && A->dy[i] == 0)
        {
          L->dx[i] = A->dx[i];
          L->dy[i] = -A->dx[i];
        }
        else if (A->dx[i] == 0 && abs(A->dy[i]) == 1)
        {
          L->dx[i] = A->dy[i];
          L->dy[i] = A->dy[i];
        }
        else if (A->dx[i] == 1 && A->dy[i] == -1)
        {
          L->dx[i] = 0;
          L->dy[i] = -A->dx[i];
        }
        else if (A->dx[i] == -1 && A->dy[i] == 1)
        {
          L->dx[i] = 0;
          L->dy[i] = -A->dx[i];
        }
        else if (A->dx[i] == 1 && A->dy[i] == 1)
        {
          L->dx[i] = A->dx[i];
          L->dy[i] = 0;
        }
        else if (A->dx[i] == -1 && A->dy[i] == -1)
        {
          L->dx[i] = A->dx[i];
          L->dy[i] = 0;
        }
      }
      return (L);
    }

  } // namespace AdjRel
} // namespace gft

namespace gft
{

  namespace Image32
  {

    sImage32 *Render(sAdjRel *A)
    {
      int i, dx, dy, dxmax, dymax;
      sImage32 *mask;
      Pixel v, u;

      dxmax = dymax = 0;
      for (i = 0; i < A->n; i++)
      {
        dx = abs(A->dx[i]);
        dy = abs(A->dy[i]);
        if (dx > dxmax)
          dxmax = dx;
        if (dy > dymax)
          dymax = dy;
      }

      mask = Create(dxmax * 2 + 1,
                    dymax * 2 + 1);
      u.x = dxmax;
      u.y = dymax;
      for (i = 0; i < A->n; i++)
      {
        v.x = u.x + A->dx[i];
        v.y = u.y + A->dy[i];

        if (IsValidPixel(mask, v.x, v.y))
          mask->array[v.y][v.x] = 1;
      }
      return mask;
    }

    void DrawAdjRel(sImage32 *img,
                    sAdjRel *A,
                    int p, int val)
    {
      Pixel u, v;
      int i;
      u.x = p % img->ncols;
      u.y = p / img->ncols;
      for (i = 0; i < A->n; i++)
      {
        v.x = u.x + A->dx[i];
        v.y = u.y + A->dy[i];
        if (IsValidPixel(img, v.x, v.y))
          img->array[v.y][v.x] = val;
      }
    }

    sImage32 *GetBoundaries(sImage32 *img,
                            sAdjRel *A)
    {
      sImage32 *himg = NULL;
      int p, q, i;
      Pixel u, v;

      himg = Create(img);

      for (u.y = 0; u.y < himg->nrows; u.y++)
      {
        for (u.x = 0; u.x < himg->ncols; u.x++)
        {
          p = u.x + u.y * himg->ncols;
          if (img->data[p] != 0)
          {
            for (i = 1; i < A->n; i++)
            {
              v.x = u.x + A->dx[i];
              v.y = u.y + A->dy[i];
              if (IsValidPixel(himg, v.x, v.y))
              {
                q = v.x + v.y * himg->ncols;
                if (img->data[p] != img->data[q])
                {
                  himg->data[p] = img->data[p];
                  break;
                }
              }
            }
          }
        }
      }
      return (himg);
    }

  } // namespace Image32

} // namespace gft

namespace gft
{
  namespace CImage
  {

    void DrawAdjRel(sCImage *cimg,
                    sAdjRel *A,
                    int p, int color)
    {
      Pixel u, v;
      int i;
      u.x = p % (cimg->C[0]->ncols);
      u.y = p / (cimg->C[0]->ncols);
      for (i = 0; i < A->n; i++)
      {
        v.x = u.x + A->dx[i];
        v.y = u.y + A->dy[i];
        if (gft::Image32::IsValidPixel(cimg->C[0], v.x, v.y))
        {
          cimg->C[0]->array[v.y][v.x] = gft::Color::Channel0(color);
          cimg->C[1]->array[v.y][v.x] = gft::Color::Channel1(color);
          cimg->C[2]->array[v.y][v.x] = gft::Color::Channel2(color);
        }
      }
    }

  } // namespace CImage
} // namespace gft
