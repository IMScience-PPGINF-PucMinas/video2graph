
#include "gft_superpixels.h"
#include <queue>

namespace gft
{
	namespace Superpixels
	{

		float ColorSquaredDistance(sCImage32f *cimg_lab,
								   const sSPixelSeed &seed, int q)
		{
			float dc;
			dc = (SQUARE(cimg_lab->C[0]->data[q] - seed.actual.l) +
				  SQUARE(cimg_lab->C[1]->data[q] - seed.actual.a) +
				  SQUARE(cimg_lab->C[2]->data[q] - seed.actual.b));
			return dc;
		}

		float LastColorSquaredDistance(const sSPixelSeed &seed)
		{
			float dc;
			dc = (SQUARE(seed.last_computed.l - seed.actual.l) +
				  SQUARE(seed.last_computed.a - seed.actual.a) +
				  SQUARE(seed.last_computed.b - seed.actual.b));
			return dc;
		}

		float LastSpatialSquaredDistance(const sSPixelSeed &seed)
		{
			float ds;
			ds = (SQUARE(seed.last_computed.i - seed.actual.i) +
				  SQUARE(seed.last_computed.j - seed.actual.j));
			return ds;
		}

		void removeSubTree(int q_in,
						   sImage32 *label,
						   sHeap *Q,
						   sAdjPxl *P,
						   sImage32 *pred,
						   float *cost,
						   sAdjRel *A)
		{
			std::queue<int> path, frontier_path;
			int i, p, q;

			path.push(q_in);
			frontier_path.push(q_in);

			while (!path.empty())
			{
				p = path.front();
				path.pop();

				label->data[p] = NIL;
				pred->data[p] = NIL;

				if (Q->color[p] == GRAY)
				{
					gft::Heap::Delete_MinPolicy(Q, p);
					Q->cost[p] = FLT_MAX;
				}
				else
				{
					Q->color[p] = WHITE;
					Q->cost[p] = FLT_MAX;
				}

				for (i = 1; i < P->n; i++)
				{
					q = p + P->dp[i];
					if (p == pred->data[q])
						path.push(q);
					else if (label->data[q] != -2)
						frontier_path.push(q);
				}
			}

			while (!frontier_path.empty())
			{
				p = frontier_path.front();
				frontier_path.pop();
				if (Q->cost[p] != FLT_MAX)
				{
					if (Q->color[p] != BLACK)
					{
						gft::Heap::Update_MinPolicy(Q, p, cost[p]);
					}
					else
					{
						gft::Heap::Insert_MinPolicy(Q, p);
					}
				}
			}
		}

		void RunSPixelsByIFT(sCImage32f *cimg_lab,
							 sImage32 *label,
							 sHeap *Q,
							 sAdjRel *A,
							 sAdjPxl *P,
							 sImage32 *pred,
							 float *cost,
							 float alpha,
							 //float beta,
							 sSPixelSeed *seeds, int nseeds)
		{
			float *dpq = NULL;
			int i, p, q, s;
			float wl, wg, tmp, alpha_2;

			alpha_2 = alpha * alpha;
			dpq = (float *)AllocFloatArray(A->n);
			for (i = 1; i < A->n; i++)
			{
				dpq[i] = sqrtf(A->dx[i] * A->dx[i] +
							   A->dy[i] * A->dy[i]);
			}

			for (s = 0; s < nseeds; s++)
			{
				p = seeds[s].actual.j + label->ncols * seeds[s].actual.i;
				label->data[p] = s;
				pred->data[p] = NIL;
				cost[p] = 0.0;
				Heap::Update_MinPolicy(Q, p, 0.0);
			}

			while (!Heap::IsEmpty(Q))
			{
				Heap::Remove_MinPolicy(Q, &p);
				s = label->data[p];

				for (i = 1; i < P->n; i++)
				{
					q = p + P->dp[i];

					if (Q->color[q] != BLACK)
					{
						wl = alpha_2 * ColorSquaredDistance(cimg_lab, seeds[s], q);
						wl *= (wl * wl);
						wl *= wl;
						tmp = cost[p] + wl + dpq[i];
						if (tmp < cost[q])
						{
							gft::Heap::Update_MinPolicy(Q, q, tmp);
							pred->data[q] = p;
							label->data[q] = s;
						}
					}
				}
			}

			for (s = 0; s < nseeds; s++)
			{
				seeds[s].last_computed = seeds[s].actual;
				seeds[s].recompute = false;
				seeds[s].updateseeds = true;
			}
			FreeFloatArray(&dpq);
		}

		void RunSPixelsByDIFT(sCImage32f *cimg_lab,
							  sImage32 *label,
							  sHeap *Q,
							  sAdjRel *A,
							  sAdjPxl *P,
							  sImage32 *pred,
							  float *cost,
							  float alpha,
							  //float beta,
							  sSPixelSeed *seeds, int nseeds)
		{
			float *dpq = NULL;
			int i, p, q, s;
			float wl, wg, tmp, alpha_2;

			alpha_2 = alpha * alpha;
			dpq = (float *)AllocFloatArray(A->n);
			for (i = 1; i < A->n; i++)
			{
				dpq[i] = sqrtf(A->dx[i] * A->dx[i] +
							   A->dy[i] * A->dy[i]);
			}

			for (s = 0; s < nseeds; s++)
			{
				if (seeds[s].recompute)
				{
					p = seeds[s].actual.j + label->ncols * seeds[s].actual.i;
					label->data[p] = s;
					pred->data[p] = NIL;
					cost[p] = 0.0;
					Heap::Update_MinPolicy(Q, p, 0.0);
					seeds[s].updateseeds = true;
				}
			}

			while (!Heap::IsEmpty(Q))
			{
				Heap::Remove_MinPolicy(Q, &p);
				s = label->data[p];

				for (i = 1; i < P->n; i++)
				{
					q = p + P->dp[i];

					if (Q->color[q] != BLACK)
					{
						wl = alpha_2 * ColorSquaredDistance(cimg_lab, seeds[s], q);
						wl *= (wl * wl);
						wl *= wl;
						tmp = cost[p] + wl + dpq[i];
						if (tmp < cost[q] || pred->data[q] == p)
						{

							if (tmp > cost[q])
							{
								seeds[label->data[q]].updateseeds = true;
								removeSubTree(q, label, Q, P, pred, cost, A);
							}
							else
							{
								if (tmp == cost[q] && pred->data[q] == p &&
									label->data[p] != label->data[q] &&
									label->data[q] != NIL)
								{
									seeds[label->data[q]].updateseeds = true;
									removeSubTree(q, label, Q, P, pred, cost, A);
								}
								else
								{
									if (label->data[q] != s)
									{
										if (label->data[q] != NIL)
											seeds[label->data[q]].updateseeds = true;
										seeds[s].updateseeds = true;
									}
									gft::Heap::Update_MinPolicy(Q, q, tmp);
									pred->data[q] = p;
									label->data[q] = s;
								}
							}
						}
					}
				}
			}

			for (s = 0; s < nseeds; s++)
			{
				if (seeds[s].recompute)
				{
					seeds[s].last_computed = seeds[s].actual;
					seeds[s].recompute = false;
				}
			}
			FreeFloatArray(&dpq);
		}

		float WeightedDistanceMeasure(sCImage32f *cimg_lab,
									  sSPixelSeed seed, int i, int j,
									  float m, float S)
		{
			float l, a, b;
			float dc, ds, D;
			l = cimg_lab->C[0]->array[i][j];
			a = cimg_lab->C[1]->array[i][j];
			b = cimg_lab->C[2]->array[i][j];
			dc = sqrtf(SQUARE(seed.actual.l - l) +
					   SQUARE(seed.actual.a - a) +
					   SQUARE(seed.actual.b - b));
			ds = sqrtf(SQUARE(seed.actual.j - j) +
					   SQUARE(seed.actual.i - i));
			D = sqrtf(dc * dc + (ds / S) * (ds / S) * m * m);
			//Na verdade, nao precisa do sqrtf pois nao vai mudar a ordem entre diferentes distancias.
			return D;
		}

		void Move2LowerGradient(sImage32 *grad,
								int *j, int *i)
		{
			sAdjRel *A;
			int u_x, u_y;
			int v_x, v_y;
			int k, g, Gmin;
			u_x = *j;
			u_y = *i;
			Gmin = grad->array[u_y][u_x];
			A = AdjRel::Circular(1.5);
			for (k = 0; k < A->n; k++)
			{
				v_x = u_x + A->dx[k];
				v_y = u_y + A->dy[k];
				if (Image32::IsValidPixel(grad, v_x, v_y))
				{
					g = grad->array[v_y][v_x];
					if (g < Gmin)
					{
						Gmin = g;
						*j = v_x;
						*i = v_y;
					}
				}
			}
			AdjRel::Destroy(&A);
		}

		void Postprocessing_CC(sImage32 *label,
							   sSPixelSeed *seeds, int nseeds)
		{
			sImage32 *CC;
			sAdjRel *A;
			int p, pp, q, n, i, lb = 0;
			int px, py, ppx, ppy, qx, qy;
			int adjlabel = 0;
			int SUPSZ;
			int *Q;
			int inic, fim;

			n = label->n;
			SUPSZ = n / nseeds;
			A = AdjRel::Neighborhood_4();
			A->dx[1] = -1;
			A->dy[1] = 0; /* left */
			A->dx[2] = 0;
			A->dy[2] = -1; /* top */
			A->dx[3] = 1;
			A->dy[3] = 0; /* right */
			A->dx[4] = 0;
			A->dy[4] = 1; /* bottom */
			CC = Image32::Create(label->ncols, label->nrows);
			Q = (int *)malloc(sizeof(int) * n);
			Image32::Set(CC, NIL);
			for (p = 0; p < n; p++)
			{
				if (CC->data[p] != NIL)
					continue;

				CC->data[p] = lb;
				px = p % label->ncols;
				py = p / label->ncols;
				for (i = 1; i < A->n; i++)
				{
					qx = px + A->dx[i];
					qy = py + A->dy[i];
					if (IsValidPixel(label, qx, qy))
					{
						q = qx + qy * label->ncols;
						if (CC->data[q] != NIL)
							adjlabel = CC->data[q];
					}
				}

				inic = fim = 0;
				Q[fim] = p;
				fim++;

				while (inic < fim)
				{
					pp = Q[inic];
					inic++;
					ppx = pp % label->ncols;
					ppy = pp / label->ncols;
					for (i = 1; i < A->n; i++)
					{
						qx = ppx + A->dx[i];
						qy = ppy + A->dy[i];
						if (IsValidPixel(label, qx, qy))
						{
							q = qx + qy * label->ncols;
							if (CC->data[q] == NIL &&
								label->data[pp] == label->data[q])
							{
								Q[fim] = q;
								fim++;
								CC->data[q] = lb;
							}
						}
					}
				}
				//----------------------------------------------------
				// If segment size is less then a limit, assign an
				// adjacent label found before, and decrement label count.
				//----------------------------------------------------
				if (fim <= SUPSZ / 4)
				{
					for (inic = 0; inic < fim; inic++)
					{
						pp = Q[inic];
						CC->data[pp] = adjlabel;
					}
					lb--;
				}
				lb++;
			}

			for (p = 0; p < n; p++)
			{
				label->data[p] = CC->data[p];
			}

			Image32::Destroy(&CC);
			AdjRel::Destroy(&A);
			free(Q);
		}

		void Postprocessing(sImage32 *label,
							sSPixelSeed *seeds, int nseeds)
		{
			sImage32 *tmp, *orphans;
			sQueue *Q;
			sAdjRel *A;
			int s, p, q, k;
			int p_x, p_y, q_x, q_y;
			int orphaned = 0;
			A = AdjRel::Neighborhood_4();
			Q = Queue::Create(label->n);
			tmp = Image32::Create(label->ncols, label->nrows);
			orphans = Image32::Create(label->ncols, label->nrows);
			Image32::Set(tmp, NIL);

			/*
	for(p = 0; p < label->n; p++){
	if(label->data[p] == NIL)
	printf("NIL\n");
	}
      */

			for (s = 0; s < nseeds; s++)
			{
				p = seeds[s].last_computed.j + label->ncols * seeds[s].last_computed.i;
				if (label->data[p] != s)
					printf("Seed %d=(%d, %d) outside its superpixel\n",
						   s, seeds[s].last_computed.j, seeds[s].last_computed.i);
				tmp->data[p] = s;
				Queue::Push(Q, p);
			}

			while (!Queue::IsEmpty(Q))
			{
				p = Queue::Pop(Q);
				p_x = p % label->ncols;
				p_y = p / label->ncols;
				for (k = 1; k < A->n; k++)
				{
					q_x = p_x + A->dx[k];
					q_y = p_y + A->dy[k];
					if (Image32::IsValidPixel(label, q_x, q_y))
					{
						q = q_x + q_y * label->ncols;
						if (label->data[p] == label->data[q] && tmp->data[q] == NIL)
						{
							tmp->data[q] = label->data[p];
							Queue::Push(Q, q);
						}
					}
				}
			}

			for (p = 0; p < label->n; p++)
			{
				if (tmp->data[p] == NIL)
				{
					orphans->data[p] = 255;
					orphaned++;
				}
			}

			Image32::Write(orphans, (char *)"orphans.pgm");
			printf("orphaned: %d\n", orphaned);

			gft::Queue::Reset(Q);

			for (p = 0; p < label->n; p++)
			{
				if (tmp->data[p] == NIL)
					continue;
				p_x = p % label->ncols;
				p_y = p / label->ncols;
				for (k = 1; k < A->n; k++)
				{
					q_x = p_x + A->dx[k];
					q_y = p_y + A->dy[k];
					if (Image32::IsValidPixel(label, q_x, q_y))
					{
						if (tmp->array[q_y][q_x] == NIL)
						{
							Queue::Push(Q, p);
							break;
						}
					}
				}
			}

			while (!Queue::IsEmpty(Q))
			{
				p = Queue::Pop(Q);
				p_x = p % label->ncols;
				p_y = p / label->ncols;
				for (k = 1; k < A->n; k++)
				{
					q_x = p_x + A->dx[k];
					q_y = p_y + A->dy[k];
					if (Image32::IsValidPixel(label, q_x, q_y))
					{
						q = q_x + q_y * label->ncols;
						if (tmp->data[q] == NIL)
						{
							tmp->data[q] = tmp->data[p];
							Queue::Push(Q, q);
						}
					}
				}
			}

			for (p = 0; p < label->n; p++)
			{
				label->data[p] = tmp->data[p];
			}

			AdjRel::Destroy(&A);
			Image32::Destroy(&tmp);
			Image32::Destroy(&orphans);
			Queue::Destroy(&Q);
		}

		sSPixelSeed *GetSPixelsSeeds(sCImage32f *cimg_lab,
									 sImage32 *grad,
									 int k,
									 int *nseeds)
		{
			int N, ncols, nrows, i, j;
			int S;
			sSPixelSeed *seeds;
			int s;
			int xstrips, ystrips;
			int xerr, yerr;
			double xerrperstrip, yerrperstrip;
			int xoff, yoff;
			int ye, xe;
			int perturbseeds = false;

			*nseeds = 0;

			ncols = cimg_lab->C[0]->ncols;
			nrows = cimg_lab->C[0]->nrows;
			N = cimg_lab->C[0]->n;
			S = ROUND(sqrtf((float)N / (float)k));

			xstrips = (0.5 + double(ncols) / double(S));
			ystrips = (0.5 + double(nrows) / double(S));

			xerr = ncols - S * xstrips;
			if (xerr < 0)
			{
				xstrips--;
				xerr = ncols - S * xstrips;
			}
			yerr = nrows - S * ystrips;
			if (yerr < 0)
			{
				ystrips--;
				yerr = nrows - S * ystrips;
			}

			xerrperstrip = double(xerr) / double(xstrips);
			yerrperstrip = double(yerr) / double(ystrips);

			//label = Image32::Create(ncols, nrows);

			xoff = S / 2;
			yoff = S / 2;

			*nseeds = xstrips * ystrips;

			seeds = (sSPixelSeed *)malloc((*nseeds) * sizeof(sSPixelSeed));

			s = 0;
			for (i = 0; i < ystrips; i++)
			{
				ye = i * yerrperstrip;
				for (j = 0; j < xstrips; j++)
				{
					xe = j * xerrperstrip;
					seeds[s].actual.j = (j * S + xoff + xe);
					seeds[s].actual.i = (i * S + yoff + ye);

					//printf("Seeds: (%d,%d)\n", seeds[s].actual.j, seeds[s].actual.i);

					//if(Image32::IsValidPixel(label, seeds[s].actual.j, seeds[s].actual.i)){
					//label->array[seeds[s].actual.i][seeds[s].actual.j] = 128;
					//Image32::DrawCircle(label, seeds[s].actual.j, seeds[s].actual.i, 4.5, 128);
					//}
					s++;
				}
			}

			for (s = 0; s < *nseeds; s++)
			{
				if (perturbseeds && grad != NULL)
					Move2LowerGradient(grad, &(seeds[s].actual.j), &(seeds[s].actual.i));
				i = seeds[s].actual.i;
				j = seeds[s].actual.j;
				//label->array[i][j] = 255;
				seeds[s].actual.l = cimg_lab->C[0]->array[i][j];
				seeds[s].actual.a = cimg_lab->C[1]->array[i][j];
				seeds[s].actual.b = cimg_lab->C[2]->array[i][j];
				seeds[s].recompute = true;
				seeds[s].updateseeds = true;
			}

			//Image32::Write(label, (char *)"seeds.pgm");
			//Image32::Destroy(&label);

			//printf("Initial number of seeds: %d\n", *nseeds);
			return seeds;
		}

		sSPixelSeed *GetSPixelsSeedsHexagon(sCImage32f *cimg_lab,
											sImage32 *grad,
											int k,
											int *nseeds)
		{
			int N, ncols, nrows, i, j;
			sSPixelSeed *seeds;
			int s;
			int perturbseeds = false;
			float dv, dh, offset;
			*nseeds = 0;
			ncols = cimg_lab->C[0]->ncols;
			nrows = cimg_lab->C[0]->nrows;
			N = cimg_lab->C[0]->n;

			dv = sqrt((2.0 * (float)N) / ((float)k * sqrt(3.0)));
			dh = dv * sqrt(3.0);
			dv = dv / 2.0;

			s = 0;
			for (i = 0; ROUND(i * dv) < nrows; i++)
			{
				if (i % 2 == 0)
					offset = dh / 2.0;
				else
					offset = 0.0;
				for (j = 0; ROUND(j * dh + offset) < ncols; j++)
					s++;
			}

			*nseeds = s;
			seeds = (sSPixelSeed *)malloc((*nseeds) * sizeof(sSPixelSeed));

			s = 0;
			for (i = 0; ROUND(i * dv) < nrows; i++)
			{
				if (i % 2 == 0)
					offset = dh / 2.0;
				else
					offset = 0.0;

				for (j = 0; ROUND(j * dh + offset) < ncols; j++)
				{
					seeds[s].actual.j = ROUND(j * dh + offset);
					seeds[s].actual.i = ROUND(i * dv);
					s++;
				}
			}

			for (s = 0; s < *nseeds; s++)
			{
				if (perturbseeds && grad != NULL)
					Move2LowerGradient(grad, &(seeds[s].actual.j), &(seeds[s].actual.i));
				i = seeds[s].actual.i;
				j = seeds[s].actual.j;
				//label->array[i][j] = 255;
				seeds[s].actual.l = cimg_lab->C[0]->array[i][j];
				seeds[s].actual.a = cimg_lab->C[1]->array[i][j];
				seeds[s].actual.b = cimg_lab->C[2]->array[i][j];
				seeds[s].recompute = true;
				seeds[s].updateseeds = true;
			}

			return seeds;
		}

		sImage32 *IFT_SLIC(sCImage *cimg, bool colored,
						   int k, float alpha, //float beta,
						   float err_dc, float err_ds, int itMax)
		{
			sCImage32f *cimg_lab, *temp;
			sImage32 *label, *flabel;
			sHeap *Q;
			sImage32 *pred;
			float *cost;
			int N, n, ncols, nrows, p, q, i;
			sSPixelSeed *seeds;
			int it, s, nseeds = 0;
			sAdjRel *A;
			sAdjPxl *P;
			//---------------------
			//char filename[512];
			clock_t end, start;
			double time_lab, time_proc;

			A = AdjRel::Neighborhood_4();

			start = clock();
			if (colored)
				cimg_lab = CImage32f::RGB2Lab(cimg);
			else
				cimg_lab = CImage32f::Clone(cimg);
			end = clock();
			time_lab = ((double)(end - start)) / CLOCKS_PER_SEC;
			printf("RGB2Lab Time: %f sec\n", time_lab);

			start = clock();

			N = cimg->C[0]->n;
			//seeds = GetSPixelsSeedsHexagon(cimg_lab, NULL, k, &nseeds);
			seeds = GetSPixelsSeeds(cimg_lab, NULL, k, &nseeds);

			temp = CImage32f::AddFrame(cimg_lab, 1, 0, 0, 0);
			CImage32f::Destroy(&cimg_lab);
			cimg_lab = temp;
			for (s = 0; s < nseeds; s++)
			{
				seeds[s].actual.i += 1;
				seeds[s].actual.j += 1;
			}

			n = cimg_lab->C[0]->n;
			ncols = cimg_lab->C[0]->ncols;
			nrows = cimg_lab->C[0]->nrows;

			P = AdjRel::AdjPixels(A, ncols);
			label = Image32::Create(ncols, nrows);
			pred = Image32::Create(ncols, nrows);
			Image32::Set(pred, NIL);
			Image32::Set(label, NIL);
			cost = AllocFloatArray(n);
			Q = Heap::Create(n, cost);
			for (p = 0; p < n; p++)
			{
				cost[p] = FLT_MAX;
			}
			//Set Frame:
			for (p = 0; p < ncols; p++)
			{
				Q->color[p] = BLACK;
				Q->color[p + ncols * (nrows - 1)] = BLACK;
				label->data[p] = -2;
				label->data[p + ncols * (nrows - 1)] = -2;
			}
			for (p = 1; p < nrows - 1; p++)
			{
				Q->color[p * ncols] = BLACK;
				Q->color[p * ncols + ncols - 1] = BLACK;
				label->data[p * ncols] = -2;
				label->data[p * ncols + ncols - 1] = -2;
			}

			for (it = 0; it < itMax; it++)
			{

				if (it == 0)
					RunSPixelsByIFT(cimg_lab, label,
									Q, A, P, pred, cost,
									alpha, //beta,
									seeds, nseeds);
				else
					RunSPixelsByDIFT(cimg_lab, label,
									 Q, A, P, pred, cost,
									 alpha, //beta,
									 seeds, nseeds);

				//---------------------
				//sprintf(filename, "label%02d.pgm", it);
				//Image32::Write(label, filename);
				//---------------------

				if (it < itMax - 1)
				{

					UpdateSPixelsSeeds(cimg_lab, label, pred, P,
									   seeds, nseeds, true);

					for (s = 0; s < nseeds; s++)
					{
						if (LastColorSquaredDistance(seeds[s]) > err_dc * err_dc ||
							LastSpatialSquaredDistance(seeds[s]) > err_ds * err_ds ||
							seeds[s].n < (N / nseeds) / 5.0)
						{
							seeds[s].recompute = true;
						}
						//#############
						else
						{
							seeds[s].actual = seeds[s].last_computed;
							seeds[s].recompute = false;
						}
						//#############
					}

					for (p = 0; p < n; p++)
					{
						s = label->data[p];
						if (s < 0)
							continue;
						Q->color[p] = WHITE;
						if (seeds[s].recompute)
						{
							cost[p] = FLT_MAX;
							pred->data[p] = NIL;
							label->data[p] = NIL;
						}
					}

					//Percorre pixels:
					//Para pixels com NIL, verifica vizinhos,
					//Se tem vizinho != NIL, entao insere na fila esse vizinho
					//com o custo que ele tinha.
					for (p = 0; p < n; p++)
					{
						if (label->data[p] != NIL)
							continue;
						for (i = 1; i < P->n; i++)
						{
							q = p + P->dp[i];
							if (label->data[q] >= 0)
								Heap::Update_MinPolicy(Q, q, cost[q]);
						}
					}
				}
			}

			Heap::Destroy(&Q);
			Image32::Destroy(&pred);
			FreeFloatArray(&cost);
			CImage32f::Destroy(&cimg_lab);
			AdjRel::Destroy(&A);
			AdjRel::DestroyAdjPxl(&P);
			free(seeds);
			flabel = Image32::RemFrame(label, 1);
			Image32::Destroy(&label);

			end = clock();
			time_proc = ((double)(end - start)) / CLOCKS_PER_SEC;
			printf("Proc Time: %f sec\n", time_proc);

			printf("IFT-SLIC Time: %f sec\n", time_proc + time_lab);

			return flabel;
		}

		sImage32 *IFT_SLIC(sCImage *cimg,
						   int k, float alpha, //float beta,
						   float err_dc, float err_ds, int itMax)
		{
			return IFT_SLIC(cimg, true,
							k, alpha, //beta,
							err_dc, err_ds, itMax);
		}

		sImage32 *IFT_SLIC(sImage32 *img,
						   int k, float alpha, //float beta,
						   float err_dc, float err_ds, int itMax)
		{
			sCImage *cimg;
			sImage32 *label;
			cimg = CImage::Clone(img);
			label = IFT_SLIC(cimg, false,
							 k, alpha / sqrtf(3.0), //beta,
							 err_dc, err_ds, itMax);
			CImage::Destroy(&cimg);
			return label;
		}

		void UpdateSPixelsSeeds(sCImage32f *cimg_lab,
								sImage32 *label,
								sImage32 *pred,
								sAdjPxl *P,
								sSPixelSeed *seeds,
								int nseeds,
								int inside)
		{
			int s, p, q, N, i, j, k, ci, cj, cd, d;
			int ncols, nrows;
			sQueue *Q;
			ncols = cimg_lab->C[0]->ncols;
			nrows = cimg_lab->C[0]->nrows;
			N = cimg_lab->C[0]->n;
			for (s = 0; s < nseeds; s++)
			{
				if (seeds[s].updateseeds)
				{
					seeds[s].n = 0;
					seeds[s].actual.i = 0;
					seeds[s].actual.j = 0;
					seeds[s].actual.l = 0.0;
					seeds[s].actual.a = 0.0;
					seeds[s].actual.b = 0.0;
				}
			}

			for (p = 0; p < N; p++)
			{
				s = label->data[p];
				if (s < 0)
					continue;
				if (seeds[s].updateseeds)
				{
					i = p / ncols;
					j = p % ncols;
					seeds[s].n++;
					seeds[s].actual.i += i;
					seeds[s].actual.j += j;
					seeds[s].actual.l += cimg_lab->C[0]->data[p];
					seeds[s].actual.a += cimg_lab->C[1]->data[p];
					seeds[s].actual.b += cimg_lab->C[2]->data[p];
				}
			}

			for (s = 0; s < nseeds; s++)
			{
				if (seeds[s].updateseeds)
				{
					if (seeds[s].n > 0)
					{
						seeds[s].actual.i /= seeds[s].n;
						seeds[s].actual.j /= seeds[s].n;
						seeds[s].actual.l /= seeds[s].n;
						seeds[s].actual.a /= seeds[s].n;
						seeds[s].actual.b /= seeds[s].n;
					}
					else
						printf("s: %d at (%d, %d)\n",
							   s, seeds[s].actual.j, seeds[s].actual.i);
				}
			}

			//------------------------------------------
			//To enforce the placement of seeds inside superpixels:
			if (inside && pred != NULL && P != NULL)
			{
				for (s = 0; s < nseeds; s++)
				{
					if (!seeds[s].updateseeds)
						continue;
					p = seeds[s].actual.i * ncols + seeds[s].actual.j;
					if (label->data[p] != s)
					{
						Q = Queue::Create(seeds[s].n);
						cj = seeds[s].last_computed.j;
						ci = seeds[s].last_computed.i;
						cd = SQUARE(cj - seeds[s].actual.j) + SQUARE(ci - seeds[s].actual.i);
						p = ci * ncols + cj;
						Queue::Push(Q, p);
						while (!Queue::IsEmpty(Q))
						{
							p = Queue::Pop(Q);
							for (k = 1; k < P->n; k++)
							{
								q = p + P->dp[k];
								if (pred->data[q] == p)
								{
									Queue::Push(Q, q);
									i = q / ncols;
									j = q % ncols;
									d = SQUARE(j - seeds[s].actual.j) + SQUARE(i - seeds[s].actual.i);
									if (d < cd)
									{
										cj = j;
										ci = i;
										cd = d;
									}
								}
							}
						}
						Queue::Destroy(&Q);
						seeds[s].actual.j = cj;
						seeds[s].actual.i = ci;
					}
				}
			}
			for (s = 0; s < nseeds; s++)
				seeds[s].updateseeds = false;
		}

		sImage32 *mySLIC(sCImage *cimg, bool colored,
						 int k, float m)
		{
			sCImage32f *cimg_lab;
			sImage32 *label, *img, *grad;
			float *dist, D;
			int N, ncols, nrows, p, i, j, si, sj;
			int S;
			sSPixelSeed *seeds;
			int it, s, nseeds = 0;

			if (colored)
			{
				cimg_lab = CImage32f::RGB2Lab(cimg);
				img = CImage::Luminosity(cimg);
			}
			else
			{
				cimg_lab = CImage32f::Clone(cimg);
				img = Image32::Clone(cimg->C[0]);
			}
			//Image32::Write(img, (char *)"luminosity.pgm");
			grad = Image32::SobelFilter(img);
			//Image32::Write(grad, (char *)"grad.pgm");

			ncols = cimg->C[0]->ncols;
			nrows = cimg->C[0]->nrows;
			N = cimg->C[0]->n;
			S = ROUND(sqrtf((float)N / (float)k));

			seeds = GetSPixelsSeeds(cimg_lab, grad, k, &nseeds);

			Image32::Destroy(&img);
			Image32::Destroy(&grad);

			label = Image32::Create(ncols, nrows);
			dist = AllocFloatArray(N);

			for (it = 0; it < 10; it++)
			{

				Image32::Set(label, NIL);
				for (p = 0; p < N; p++)
					dist[p] = FLT_MAX;

				for (s = 0; s < nseeds; s++)
				{
					si = seeds[s].actual.i;
					sj = seeds[s].actual.j;
					p = sj + si * ncols;
					dist[p] = 0;
					label->data[p] = s;

					for (i = ROUND(si - S); i <= ROUND(si + S); i++)
					{
						for (j = ROUND(sj - S); j <= ROUND(sj + S); j++)
						{
							if (Image32::IsValidPixel(label, j, i))
							{
								D = WeightedDistanceMeasure(cimg_lab, seeds[s], i, j, m, S);
								p = j + i * ncols;
								if (D < dist[p])
								{
									dist[p] = D;
									label->data[p] = s;
								}
							}
						}
					}
				}

				UpdateSPixelsSeeds(cimg_lab, label, NULL, NULL,
								   seeds, nseeds, false);
			}

			//Image32::Write(label, (char *)"label_before.pgm");
			Postprocessing_CC(label, seeds, nseeds);

			CImage32f::Destroy(&cimg_lab);
			FreeFloatArray(&dist);
			free(seeds);
			return label;
		}

		sImage32 *mySLIC(sCImage *cimg,
						 int k, float m)
		{
			return mySLIC(cimg, true, k, m);
		}

		sImage32 *mySLIC(sImage32 *img,
						 int k, float m)
		{
			sCImage *cimg;
			sImage32 *label;
			cimg = CImage::Clone(img);
			label = mySLIC(cimg, false, k, m * sqrtf(3.0));
			CImage::Destroy(&cimg);
			return label;
		}

		int GetNumberOfSuperPixels(sImage32 *label)
		{
			int p;
			int Lmin, Lmax;
			Lmin = INT_MAX;
			Lmax = INT_MIN;
			for (p = 0; p < label->n; p++)
			{
				if (label->data[p] < Lmin)
					Lmin = label->data[p];
				if (label->data[p] > Lmax)
					Lmax = label->data[p];
			}
			//printf("Lmin: %d, Lmax: %d\n", Lmin, Lmax);
			return Lmax - Lmin + 1;
		}

		//-----------------------------------------------

		float FeatureSquaredDistance(sScene32 *scn,
									 const sSVoxelSeed &seed, int q)
		{
			float df;
			df = SQUARE(scn->data[q] - seed.actual.l);
			//df = fabsf(scn->data[q] - seed.actual.l);
			return df;
		}

		float LastFeatureSquaredDistance(const sSVoxelSeed &seed)
		{
			float dc;
			dc = SQUARE(seed.last_computed.l - seed.actual.l);
			//dc = fabsf(seed.last_computed.l - seed.actual.l);
			return dc;
		}

		float LastSpatialSquaredDistance(const sSVoxelSeed &seed)
		{
			float ds;
			ds = (SQUARE(seed.last_computed.i - seed.actual.i) +
				  SQUARE(seed.last_computed.j - seed.actual.j) +
				  SQUARE(seed.last_computed.k - seed.actual.k));
			return ds;
		}

		void removeSubTree3(int q_in,
							sScene32 *label,
							sHeap *Q,
							sAdjVxl *P,
							sScene32 *pred,
							float *cost,
							sAdjRel3 *A)
		{
			std::queue<int> path, frontier_path;
			int i, p, q;

			path.push(q_in);
			frontier_path.push(q_in);

			while (!path.empty())
			{
				p = path.front();
				path.pop();

				label->data[p] = NIL;
				pred->data[p] = NIL;

				if (Q->color[p] == GRAY)
				{
					gft::Heap::Delete_MinPolicy(Q, p);
					Q->cost[p] = FLT_MAX;
				}
				else
				{
					Q->color[p] = WHITE;
					Q->cost[p] = FLT_MAX;
				}

				for (i = 1; i < P->n; i++)
				{
					q = p + P->dp[i];
					if (p == pred->data[q])
						path.push(q);
					else if (label->data[q] != -2)
						frontier_path.push(q);
				}
			}

			while (!frontier_path.empty())
			{
				p = frontier_path.front();
				frontier_path.pop();
				if (Q->cost[p] != FLT_MAX)
				{
					if (Q->color[p] != BLACK)
					{
						gft::Heap::Update_MinPolicy(Q, p, cost[p]);
					}
					else
					{
						gft::Heap::Insert_MinPolicy(Q, p);
					}
				}
			}
		}

		void RunSVoxelsByIFT(sScene32 *scn,
							 sScene32 *label,
							 sHeap *Q,
							 sAdjRel3 *A,
							 sAdjVxl *V,
							 sScene32 *pred,
							 float *cost,
							 float alpha,
							 //float beta,
							 sSVoxelSeed *seeds, int nseeds)
		{
			float *dpq = NULL;
			int i, p, q, s;
			float wl, tmp, alpha_2;

			alpha_2 = alpha * alpha;
			dpq = (float *)AllocFloatArray(A->n);
			for (i = 1; i < A->n; i++)
			{
				dpq[i] = sqrtf(A->d[i].axis.x * A->d[i].axis.x +
							   A->d[i].axis.y * A->d[i].axis.y +
							   A->d[i].axis.z * A->d[i].axis.z);
			}

			for (s = 0; s < nseeds; s++)
			{
				p = Scene32::GetVoxelAddress(scn,
											 seeds[s].actual.j,
											 seeds[s].actual.i,
											 seeds[s].actual.k);
				label->data[p] = s;
				pred->data[p] = NIL;
				cost[p] = 0.0;
				Heap::Update_MinPolicy(Q, p, 0.0);
			}

			while (!Heap::IsEmpty(Q))
			{
				Heap::Remove_MinPolicy(Q, &p);
				s = label->data[p];

				for (i = 1; i < V->n; i++)
				{
					q = p + V->dp[i];

					if (Q->color[q] != BLACK)
					{
						wl = alpha_2 * FeatureSquaredDistance(scn, seeds[s], q);
						wl *= (wl * wl);
						wl *= wl;
						tmp = cost[p] + wl + dpq[i];
						if (tmp < cost[q])
						{
							Heap::Update_MinPolicy(Q, q, tmp);
							pred->data[q] = p;
							label->data[q] = s;
						}
					}
				}
			}

			for (s = 0; s < nseeds; s++)
			{
				seeds[s].last_computed = seeds[s].actual;
				seeds[s].recompute = false;
				seeds[s].updateseeds = true;
			}
			FreeFloatArray(&dpq);
		}

		void RunSVoxelsByDIFT(sScene32 *scn,
							  sScene32 *label,
							  sHeap *Q,
							  sAdjRel3 *A,
							  sAdjVxl *V,
							  sScene32 *pred,
							  float *cost,
							  float alpha,
							  //float beta,
							  sSVoxelSeed *seeds, int nseeds)
		{
			float *dpq = NULL;
			int i, p, q, s;
			float wl, tmp, alpha_2;

			alpha_2 = alpha * alpha;
			dpq = (float *)AllocFloatArray(A->n);
			for (i = 1; i < A->n; i++)
			{
				dpq[i] = sqrtf(A->d[i].axis.x * A->d[i].axis.x +
							   A->d[i].axis.y * A->d[i].axis.y +
							   A->d[i].axis.z * A->d[i].axis.z);
			}

			for (s = 0; s < nseeds; s++)
			{
				if (seeds[s].recompute)
				{
					p = Scene32::GetVoxelAddress(scn,
												 seeds[s].actual.j,
												 seeds[s].actual.i,
												 seeds[s].actual.k);
					label->data[p] = s;
					pred->data[p] = NIL;
					cost[p] = 0.0;
					Heap::Update_MinPolicy(Q, p, 0.0);
					seeds[s].updateseeds = true;
				}
			}

			while (!Heap::IsEmpty(Q))
			{
				Heap::Remove_MinPolicy(Q, &p);
				s = label->data[p];

				for (i = 1; i < V->n; i++)
				{
					q = p + V->dp[i];

					if (Q->color[q] != BLACK)
					{
						wl = alpha_2 * FeatureSquaredDistance(scn, seeds[s], q);
						wl *= (wl * wl);
						wl *= wl;
						tmp = cost[p] + wl + dpq[i];
						if (tmp < cost[q] || pred->data[q] == p)
						{

							if (tmp > cost[q])
							{
								seeds[label->data[q]].updateseeds = true;
								removeSubTree3(q, label, Q, V, pred, cost, A);
							}
							else
							{
								if (tmp == cost[q] && pred->data[q] == p &&
									label->data[p] != label->data[q] &&
									label->data[q] != NIL)
								{
									seeds[label->data[q]].updateseeds = true;
									removeSubTree3(q, label, Q, V, pred, cost, A);
								}
								else
								{
									if (label->data[q] != s)
									{
										if (label->data[q] != NIL)
											seeds[label->data[q]].updateseeds = true;
										seeds[s].updateseeds = true;
									}
									gft::Heap::Update_MinPolicy(Q, q, tmp);
									pred->data[q] = p;
									label->data[q] = s;
								}
							}
						}
					}
				}
			}

			for (s = 0; s < nseeds; s++)
			{
				if (seeds[s].recompute)
				{
					seeds[s].last_computed = seeds[s].actual;
					seeds[s].recompute = false;
				}
			}
			FreeFloatArray(&dpq);
		}

		sScene32 *IFT_SLIC(sScene32 *scn,
						   int k, float alpha, //float beta,
						   float err_df, float err_ds, int itMax)
		{
			sScene32 *label, *flabel, *pred, *fscn;
			sHeap *Q;
			float *cost;
			int N, n, xsize, ysize, zsize, p, q, i;
			sSVoxelSeed *seeds;
			int x, y, z, it, s, nseeds = 0;
			sAdjRel3 *A;
			sAdjVxl *V;
			//------------------
			//FILE *fp;

			A = AdjRel3::Spheric(1.0);

			N = scn->n;
			seeds = GetSVoxelsSeeds(scn, k, &nseeds);

			printf("Numero seeds IFT: %d\n", nseeds);
			fscn = Scene32::AddFrame(scn, 1, 0);
			for (s = 0; s < nseeds; s++)
			{
				seeds[s].actual.i += 1;
				seeds[s].actual.j += 1;
				seeds[s].actual.k += 1;
			}

			V = AdjRel3::AdjVoxels(A, fscn);

			n = fscn->n;
			xsize = fscn->xsize;
			ysize = fscn->ysize;
			zsize = fscn->zsize;
			label = Scene32::Create(fscn);
			pred = Scene32::Create(fscn);
			Scene32::Fill(pred, NIL);
			Scene32::Fill(label, NIL);
			cost = AllocFloatArray(n);
			Q = Heap::Create(n, cost);
			for (p = 0; p < n; p++)
			{
				cost[p] = FLT_MAX;
			}
			//Set Frame:
			for (p = 0; p < xsize * ysize; p++)
			{
				Q->color[p] = BLACK;
				Q->color[p + (xsize * ysize) * (zsize - 1)] = BLACK;
				label->data[p] = -2;
				label->data[p + (xsize * ysize) * (zsize - 1)] = -2;
			}
			for (z = 0; z < zsize; z++)
			{
				for (x = 0; x < xsize; x++)
				{
					y = 0;
					p = gft::Scene32::GetVoxelAddress(fscn, x, y, z);
					Q->color[p] = BLACK;
					label->data[p] = -2;
					y = ysize - 1;
					p = gft::Scene32::GetVoxelAddress(fscn, x, y, z);
					Q->color[p] = BLACK;
					label->data[p] = -2;
				}
			}
			printf("Numero zsize IFT: %d\n", zsize);
			for (z = 0; z < zsize; z++)
			{
				for (y = 0; y < ysize; y++)
				{
					x = 0;
					p = gft::Scene32::GetVoxelAddress(fscn, x, y, z);
					Q->color[p] = BLACK;
					label->data[p] = -2;
					x = xsize - 1;
					p = gft::Scene32::GetVoxelAddress(fscn, x, y, z);
					Q->color[p] = BLACK;
					label->data[p] = -2;
				}
			}

			for (it = 0; it < itMax; it++)
			{

				if (it == 0)
					RunSVoxelsByIFT(fscn, label,
									Q, A, V, pred, cost,
									alpha, //beta,
									seeds, nseeds);
				else
					RunSVoxelsByDIFT(fscn, label,
									 Q, A, V, pred, cost,
									 alpha, //beta,
									 seeds, nseeds);

				if (it < itMax - 1)
				{

					UpdateSVoxelsSeeds(fscn, label, pred, V,
									   seeds, nseeds, true);

					for (s = 0; s < nseeds; s++)
					{
						if (LastFeatureSquaredDistance(seeds[s]) > err_df * err_df ||
							LastSpatialSquaredDistance(seeds[s]) > err_ds * err_ds ||
							seeds[s].n < (N / nseeds) / 5.0)
						{
							seeds[s].recompute = true;
						}
						//#############
						else
						{
							seeds[s].actual = seeds[s].last_computed;
							seeds[s].recompute = false;
						}
						//#############
					}

					for (p = 0; p < n; p++)
					{
						s = label->data[p];
						if (s < 0)
							continue;
						Q->color[p] = WHITE;
						if (seeds[s].recompute)
						{
							cost[p] = FLT_MAX;
							pred->data[p] = NIL;
							label->data[p] = NIL;
						}
					}

					//Percorre pixels:
					//Para pixels com NIL, verifica vizinhos,
					//Se tem vizinho != NIL, entao insere na fila esse vizinho
					//com o custo que ele tinha.
					for (p = 0; p < n; p++)
					{
						if (label->data[p] != NIL)
							continue;
						for (i = 1; i < V->n; i++)
						{
							q = p + V->dp[i];
							if (label->data[q] >= 0)
								Heap::Update_MinPolicy(Q, q, cost[q]);
						}
					}
				}
			}

			Heap::Destroy(&Q);
			Scene32::Destroy(&pred);
			FreeFloatArray(&cost);
			AdjRel3::Destroy(&A);
			AdjRel3::DestroyAdjVxl(&V);
			//--------------------
			/*
      fp = fopen("SP_seeds_report.txt", "w");
      if(fp != NULL){
	fprintf(fp, "%d\n", nseeds);
	for(i = 0; i < nseeds; i++){
	  fprintf(fp, "%d %d %d %d %f\n",
		  seeds[i].actual.j,
		  seeds[i].actual.i,
		  seeds[i].actual.k,
		  i,
		  seeds[i].actual.l);
	}
	fclose(fp);
      }
      */
			//--------------------
			free(seeds);
			flabel = Scene32::RemFrame(label, 1);
			Scene32::Destroy(&label);
			return flabel;
		}

		sSVoxelSeed *GetSVoxelsSeeds(sScene32 *scn,
									 int k,
									 int *nseeds)
		{
			//Scene32::Scene32 *label;
			int N, i, j, l, S;
			sSVoxelSeed *seeds;
			int s;
			int xstrips, ystrips, zstrips;
			int xerr, yerr, zerr;
			double xerrperstrip, yerrperstrip, zerrperstrip;
			int xoff, yoff, zoff;
			int ye, xe, ze;

			*nseeds = 0;

			N = scn->n;
			if (scn->zsize > 1)
				S = ROUND(pow((double)N / (double)k, 1.0 / (1.0 + scn->zsize ))); // tava 3.0
			else
				S = ROUND(pow((double)N / (double)k, 1.0 / 2.0));
			printf("S: %d\n",S );	
			printf("scn->xsize %d, scn->ysize %d, scn->zsize: %d\n", scn->xsize, scn->ysize, scn->zsize);
			xstrips = (0.5 + double(scn->xsize) / double(S));
			ystrips = (0.5 + double(scn->ysize) / double(S));
			zstrips = (0.5 + double(scn->zsize) / double(S));
			printf("xstrips %d, ystrips %d, zstrips: %d\n", xstrips, ystrips, zstrips);
			xerr = scn->xsize - S * xstrips;
			if (xerr < 0)
			{
				xstrips--;
				xerr = scn->xsize - S * xstrips;
			}

			yerr = scn->ysize - S * ystrips;
			if (yerr < 0)
			{
				ystrips--;
				yerr = scn->ysize - S * ystrips;
			}

			zerr = scn->zsize - S * zstrips;
			if (zerr < 0)
			{
				zstrips--;
				zerr = scn->zsize - S * zstrips;
			}

			xerrperstrip = double(xerr) / double(xstrips);
			yerrperstrip = double(yerr) / double(ystrips);
			if (zstrips != 0.0)
				zerrperstrip = double(zerr) / double(zstrips);
			else
				zerrperstrip = zerr;

			//label = Scene32::Create(scn);
			printf("xstrips %d, ystrips %d, zstrips: %d  222\n", xstrips, ystrips, zstrips);
			xoff = S / 2;
			yoff = S / 2;
			if (scn->zsize > 1)
				zoff = S / 2;
			else
				zoff = 0;

			if (scn->zsize > 1)
				*nseeds = xstrips * ystrips * zstrips;
			else
			{
				zstrips = 1;
				*nseeds = xstrips * ystrips * zstrips;
			}

			seeds = (sSVoxelSeed *)malloc((*nseeds) * sizeof(sSVoxelSeed));

			s = 0;
			for (l = 0; l < zstrips; l++)
			{
				ze = l * zerrperstrip;
				for (i = 0; i < ystrips; i++)
				{
					ye = i * yerrperstrip;
					for (j = 0; j < xstrips; j++)
					{
						xe = j * xerrperstrip;

						seeds[s].actual.j = (j * S + xoff + xe);
						seeds[s].actual.i = (i * S + yoff + ye);
						seeds[s].actual.k = (l * S + zoff + ze);

						/*
	    if(Scene32::IsValidVoxel(label,
				     seeds[s].actual.j,
				     seeds[s].actual.i,
				     seeds[s].actual.k)){
	      label->array[seeds[s].actual.k][seeds[s].actual.i][seeds[s].actual.j] = 128;
	      //Image32::DrawCircle(label, seeds[s].actual.j, seeds[s].actual.i, 4.5, 128);
	    }
	    */

						s++;
					}
				}
			}

			for (s = 0; s < *nseeds; s++)
			{
				i = seeds[s].actual.i;
				j = seeds[s].actual.j;
				l = seeds[s].actual.k;
				seeds[s].actual.l = scn->array[l][i][j];
				seeds[s].recompute = true;
				seeds[s].updateseeds = true;
			}

			/*
      Scene32::Write(label, (char *)"seeds.scn");
      Scene32::Destroy(&label);
      printf("Initial number of seeds: %d\n", *nseeds);
      */
			return seeds;
		}

		void UpdateSVoxelsSeeds(sScene32 *scn,
								sScene32 *label,
								sScene32 *pred,
								sAdjVxl *V,
								sSVoxelSeed *seeds, int nseeds,
								int inside)
		{
			int s, p, q, N, i, j, k, ci, cj, ck, cd, d;
			sQueue *Q;
			N = scn->n;
			for (s = 0; s < nseeds; s++)
			{
				if (seeds[s].updateseeds)
				{
					seeds[s].n = 0;
					seeds[s].actual.i = 0;
					seeds[s].actual.j = 0;
					seeds[s].actual.k = 0;
					seeds[s].actual.l = 0.0;
				}
			}

			for (p = 0; p < N; p++)
			{
				s = label->data[p];
				if (s < 0)
					continue;
				if (seeds[s].updateseeds)
				{
					j = Scene32::GetAddressX(scn, p);
					i = Scene32::GetAddressY(scn, p);
					k = Scene32::GetAddressZ(scn, p);
					seeds[s].n++;
					seeds[s].actual.i += i;
					seeds[s].actual.j += j;
					seeds[s].actual.k += k;
					seeds[s].actual.l += scn->data[p];
				}
			}

			for (s = 0; s < nseeds; s++)
			{
				if (seeds[s].updateseeds)
				{
					if (seeds[s].n > 0)
					{
						seeds[s].actual.i /= seeds[s].n;
						seeds[s].actual.j /= seeds[s].n;
						seeds[s].actual.k /= seeds[s].n;
						seeds[s].actual.l /= seeds[s].n;
					}
					else
						printf("s: %d at (%d, %d, %d)\n",
							   s, seeds[s].actual.j, seeds[s].actual.i, seeds[s].actual.k);
				}
			}

			//------------------------------------------
			//To enforce the placement of seeds inside superpixels:
			if (inside && pred != NULL && V != NULL)
			{
				for (s = 0; s < nseeds; s++)
				{
					if (!seeds[s].updateseeds)
						continue;
					p = gft::Scene32::GetVoxelAddress(scn,
													  seeds[s].actual.j,
													  seeds[s].actual.i,
													  seeds[s].actual.k);
					if (label->data[p] != s)
					{
						Q = Queue::Create(seeds[s].n);
						cj = seeds[s].last_computed.j;
						ci = seeds[s].last_computed.i;
						ck = seeds[s].last_computed.k;
						cd = (SQUARE(cj - seeds[s].actual.j) +
							  SQUARE(ci - seeds[s].actual.i) +
							  SQUARE(ck - seeds[s].actual.k));
						p = gft::Scene32::GetVoxelAddress(scn, cj, ci, ck);
						Queue::Push(Q, p);
						while (!Queue::IsEmpty(Q))
						{
							p = Queue::Pop(Q);
							for (k = 1; k < V->n; k++)
							{
								q = p + V->dp[k];
								if (pred->data[q] == p)
								{
									Queue::Push(Q, q);
									j = Scene32::GetAddressX(scn, q);
									i = Scene32::GetAddressY(scn, q);
									k = Scene32::GetAddressZ(scn, q);
									d = (SQUARE(j - seeds[s].actual.j) +
										 SQUARE(i - seeds[s].actual.i) +
										 SQUARE(k - seeds[s].actual.k));
									if (d < cd)
									{
										cj = j;
										ci = i;
										ck = k;
										cd = d;
									}
								}
							}
						}
						Queue::Destroy(&Q);
						seeds[s].actual.j = cj;
						seeds[s].actual.i = ci;
						seeds[s].actual.k = ck;
					}
				}
			}
			for (s = 0; s < nseeds; s++)
				seeds[s].updateseeds = false;
		}

	} // namespace Superpixels
} // namespace gft
