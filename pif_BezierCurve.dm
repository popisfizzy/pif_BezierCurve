#if !defined(PIF_NOPREFIX_GENERAL) && !defined(PIF_NOPREFIX_BEZIERCURVE)
pif_BezierCurve
#else
BezierCurve
#endif
	var
		// The matrix of the current transformation. If null, it will be treated like
		// the 3x3 identity matrix.
		matrix/transform

		list
			// A list of points for the Bézier Curve. It is even in length.
			Points = new
			TransformPointCache // The points after they have been transformed by a call to
							    // transform. If this value is not null, then it's values are
							    // what will be used to get new points.

			// A cache of binomial coefficients. This is initialized during calls to methods
			// that use the _Choose() method and then cleared out afterwards.
			BinomialCache

			// A cache of the most recent call to Curve(). This will be used on subsequent
			// calls provided that (1) smoothness does not change, (2) the points do not
			// change, (3) the points are not given a new transformation, and (4) a custom
			// smoothness argument is not passed to Curve().
			CurveCache

		// The smoothness of the curve. This is how many equidistant points (when measuring
		// along the arclength) are found.
		smoothness

	New(...)
		// The argument format is x0,y0 , x1,x1, ..., xN,yN.

		if(args.len == 0)
#if !defined(PIF_NOPREFIX_GENERAL) && !defined(PIF_NOPREFIX_BEZIERCURVE)
			throw new /pif_BezierCurve/InvalidConstructorArgumentFormatException(__FILE__, __LINE__)
#else
			throw new /BezierCurve/InvalidConstructorArgumentFormatException(__FILE__, __LINE__)
#endif


		var/list/Data

		if(istype(args[1], /list))
			Data = args[1]
		else
			Data = args

		if(Data.len < 4)
#if !defined(PIF_NOPREFIX_GENERAL) && !defined(PIF_NOPREFIX_BEZIERCURVE)
			throw new /pif_BezierCurve/InvalidArgumentDataTooFewPointsException(__FILE__, __LINE__)
#else
			throw new /BezierCurve/InvalidArgumentDataTooFewPointsException(__FILE__, __LINE__)
#endif
		if((Data.len % 2) != 0)
#if !defined(PIF_NOPREFIX_GENERAL) && !defined(PIF_NOPREFIX_BEZIERCURVE)
			throw new /pif_BezierCurve/InvalidArgumentDataOddNumberException(__FILE__, __LINE__)
#else
			throw new /BezierCurve/InvalidArgumentDataOddNumberException(__FILE__, __LINE__)
#endif

		Points = Data.Copy()

	proc
		_Choose(n,k)
			// Insight for this method found at http://math.stackexchange.com/a/202559/325737

			if(BinomialCache[n] == null)
				var/list/L = new
				L.len = Length()

				BinomialCache[n] = L
			else if(BinomialCache[n][k+1] != null)
				// Because k == 0 pretty regularly, we just shift everything one to the
				// right.
				return BinomialCache[n][k+1]

			if(k > (n-k))
				k = n-k

			. = 1 // Result.
			for(var/i = 0, i < k, i ++)
				. *= (n-i)/(k-i)

			BinomialCache[n][k+1] = .
			return .

		Length()
			// Number of points, not length of the Points list, hence the division by 2.
			return Points.len / 2

		Smoothness(s_ = src.smoothness)
			if(s_ != src.smoothness)
				// We'll have to recompute the curve if the value has changed, so clear out the
				// cache.
				CurveCache = null

			smoothness = s_
			return smoothness

		Curve(_smoothness = src.smoothness)
			// Generate a list of points that we will draw lines between, in order to plot the
			// (approximated) Bezier curve.

			if(!CurveCache || !_smoothness || (_smoothness != src.smoothness))

				var
					// This is where we'll be drawing the points from.
					list/Data = TransformPointCache ? TransformPointCache : Points

					// If smoothness is not set, we'll populate the list with as many points in as
					// there are points out.
					smoothness = (_smoothness || Smoothness() || Length()) - 1

				// The output will be a list with length 2*Smoothness(), where the even indices
				// are the x coordinates and the odd indices are the y coordinates.
				CurveCache = list(
					// This is because a Bézier curve always starts on its first control point.
					// That is, if P_i is the ith control point, B(0) = P_0.
					Data[1],
					Data[2]
				)

				CurveCache.len = 2 * (smoothness + 1)

				BinomialCache = new
				BinomialCache.len = Length()

				var/N = Length()
				for(var/i = 1, i < smoothness, i ++)
					var
						delta = i / smoothness

						q = 1
						r = (1 - delta)**(N-1)

						x = 0
						y = 0

					for(var/k = 0, k < N, k ++)
						var
							px = Data[2*k+1]
							py = Data[2*k+2]

							binomial = _Choose(N-1, k)
							coefficient = q*r * binomial

						x += px * coefficient
						y += py * coefficient

						q *= delta
						r /= (1-delta)

					CurveCache[2*i+1] = x
					CurveCache[2*i+2] = y

				// A Bézier curve always ends on its last control point.
				CurveCache[CurveCache.len-1] = Data[Data.len-1]
				CurveCache[CurveCache.len  ] = Data[Data.len  ]

			return CurveCache

		ClearTransform()
			// This clears out any current transform.

			transform = null
			TransformPointCache = null
			CurveCache = null

			return 1

		Transform(matrix/M)
			// Take either a matrix, a list, or a list of arguments and applies the transformation
			// data in this object in the appropriate way.

			if(istype(M, /list))
				if(M:len == 6)
					M = new(
						M[1], M[2], M[3],
						M[4], M[5], M[6]
					)

				else
#if !defined(PIF_NOPREFIX_GENERAL) && !defined(PIF_NOPREFIX_BEZIERCURVE)
					throw new /pif_BezierCurve/InvalidArgumentFormatException(__FILE__, __LINE__)
#else
					throw new /BezierCurve/InvalidArgumentFormatException(__FILE__, __LINE__)
#endif

			else if(args.len == 6)
				M = new(
					args[1], args[2], args[3],
					args[4], args[5], args[6]
				)

			else if(args.len == 0)
				// Treat it as the identity transform.
				M = new

			else if(M == null)
				// Same.
				M = new

			else if(!istype(M))
				// Otherwise, it's an unknown argument format.
#if !defined(PIF_NOPREFIX_GENERAL) && !defined(PIF_NOPREFIX_BEZIERCURVE)
				throw new /pif_BezierCurve/InvalidArgumentFormatException(__FILE__, __LINE__)
#else
				throw new /BezierCurve/InvalidArgumentFormatException(__FILE__, __LINE__)
#endif

			// If transform is not set, we set it to M. Otherwise, we right-multiply transform by
			// M to get the new matrix.
			transform = transform ? M : (transform * M)

			// And now we re-compute the TransformPointCache list via _ComputeTransform() and
			// return its value.

			return _ComputeTransform()

		_ComputeTransform()
			// Recomputes the transform when needed, such as after applying a new transform or
			// updating the list of points.

			var
				a = transform.a
				b = transform.b
				c = transform.c
				d = transform.d
				e = transform.e
				f = transform.f

			TransformPointCache = new
			TransformPointCache.len = Points.len

			for(var/i = 1, i <= Points.len, i += 2)
				var
					x = Points[i  ]
					y = Points[i+1]

				TransformPointCache[i  ] = a*x + b*y + c
				TransformPointCache[i+1] = d*x + e*y + f

			// Since the points have changed, clear out CurveCache. This will have us re-compute
			// then curve next time Curve() is called, instead of just outputting the cached value.
			CurveCache = null

			return 1

		InsertPoint(p, x,y)
			// Insert the point (x,y) at position p. Negative arguments indicate to start adding
			// at the end.

			if(p > 0)
				Points.Insert(2*p-1, x,y)
			else
				Points.Insert(Points.len + 2*(p+1) + 1, x,y)

			_ComputeTransform()
			return Length()

		AddPoint(x,y)
			// Adds a new point to the end.
			. = InsertPoint(-1, x,y)
			_ComputeTransform()

		RemovePoint(p)
			// Removes the point at position p. Negative arguments indicate distance from the
			// end of the list.

			if(p > 0)
				Points.Cut(2*p-1, 2*p+1)
			else
				Points.Cut(Points.len + 2*(p+1) - 1, Points.len + 2*(p+1) + 1)

			_ComputeTransform()
			return Length()

		NudgePoint(p, dx,dy)
			// If (x,y) is at p, then this will do (x,y) += (dx,dy). As with InsertPoint(), when
			// p is positive it measures from the beginning of the list, and when p is negative it
			// measures from the end of the list.

			if(p > 0)
				Points[2*p-1] += dx
				Points[2*p]   += dy

			else
				Points[Points.len + 2*(p+1) - 1] += dx
				Points[Points.len + 2*(p+1)    ] += dy

			_ComputeTransform()
			return 1

		Point(p)
			// Gets the pth point. As always, positive measures from the front of the list and
			// negative measures from the back of the list.

			. = list(0, 0)

			if(p > 0)
				.[1] = Points[2*p-1]
				.[2] = Points[2*p  ]

			else
				.[1] = Points[Points.len + 2*(p+1) - 1]
				.[2] = Points[Points.len + 2*(p+1)    ]

		FirstPoint()
			return Point(1)

		LastPoint()
			return Point(-1)

		Close()
			// This will "close" the curve by adding a new point equal to the initial point,
			// provided the curve is not already closed.

			var/list
				F = FirstPoint()
				L = LastPoint()

			if( (F[1] == L[1]) && (F[2] == L[2]) )
				return 0
			else
				AddPoint(F[1], F[2])
				_ComputeTransform()
				return 1

		Open()
			// This will "open" the curve by removing the last point, but only if it's equal
			// to the initial point.

			var/list
				F = FirstPoint()
				L = LastPoint()

			if( (F[1] == L[1]) && (F[2] == L[2]) )
				RemovePoint(-1)
				_ComputeTransform()
				return 1
			else
				return 0

		ArcLength(s_ = smoothness)
			// Gives the approximate Euclidean length along the given Bézier curve, limited
			// largely by the specified smoothness. This applies to the *transformed* Bézier
			// curve if it is present, or otherwise the original.

			var/list/Curve = Curve(s_)

			. = 0
			for(var/i = 1, i < (Curve.len-1), i += 2)
				var
					dx = Curve[i+2]-Curve[i  ]
					dy = Curve[i+3]-Curve[i+1]

				. += sqrt(dx*dx + dy*dy)

		Displacement()
			// Gives the displacement from the first point to the last point, either on the transformed
			// points if they exist or on the original otherwise. This is a vector quantity.

			var
				xf = TransformPointCache ? TransformPointCache[TransformPointCache.len-1] : Points[Points.len-1]
				yf = TransformPointCache ? TransformPointCache[TransformPointCache.len  ] : Points[Points.len  ]

				xi = TransformPointCache ? TransformPointCache[1] : Points[1]
				yi = TransformPointCache ? TransformPointCache[2] : Points[2]

			return list(
				xf-xi, yf-yi
			)

		Distance()
			// Gives the distance from the first point to the last point, either on the transformed points
			// if they exist or on the original otherwise. This is a scalar quantity.

			. = Displacement()
			return sqrt(.[1]*.[1] + .[2]*.[2])

		Bezier(t)
			// This gives the point at the Bézier curve at a specific time t, rather than an approximation
			// of the curve over the range t in [0,1].

			. = list(0,0)

			var/list/Data = TransformPointCache ? TransformPointCache : Points

			if(t == 0)
				// Bezier(0) is always equal to the first point
				. = list(
					Data[1], Data[2]
				)
			else if(t == 1)
				// And Bezier(1) is always equal to the last point.
				. = list(
					Data[Data.len-1], Data[Data.len]
				)

			else
				// If it's not the first or last point, then compute it.

				var
					N = Length()

					q = 1
					r = (1 - t)**(N-1)

				BinomialCache = new
				BinomialCache.len = N

				for(var/k = 0, k < N, k ++)
					var
						px = Data[2*k+1]
						py = Data[2*k+2]

						binomial = _Choose(N-1, k)
						coefficient = q*r * binomial

					.[1] += px * coefficient
					.[2] += py * coefficient

					q *= t
					r /= (1-t)

				BinomialCache = null

		Derivative(t)
			// This gives the derivative of the Bézier curve at time t. This could be used to,
			// for example, stop an animation along a Bézier curve and move in a straight line
			// from there.

			// This code uses the following identity for d/dt Bézier(t):
			//    Bézier'(t) = n * Sum( Choose(n-1, k) * t**k * (1-t)**(n-k-1) (P(k+1) - P(k))) from k = 0 to n-1.

			. = list(0,0)

			var
				list/Data = TransformPointCache ? TransformPointCache : Points
				N = Length()-1

			if(t == 0)
				// Bezier'(0) can be simplified down to the following:
				. = list(
					N * (Data[3] - Data[1]),
					N * (Data[4] - Data[2])
				)
			else if(t == 1)
				// Bezier'(1) can be simplified down to the following:
				. = list(
					N * (Data[Data.len-1] - Data[Data.len-3]),
					N * (Data[Data.len  ] - Data[Data.len-2])
				)

			else
				// For the rest, we'll just compute the derivative directly from the previous
				// formula.

				var
					q = 1
					r = (1-t)**(N-1)

					old_x = Data[1]
					old_y = Data[2]

				BinomialCache = new
				BinomialCache.len = N

				for(var/k = 0, k < N, k ++)
					var
						new_x = Data[2*k+3]
						new_y = Data[2*k+4]

						binomial = _Choose(N-1, k)
						coefficient = N * q*r * binomial

					.[1] += coefficient * (new_x - old_x)
					.[2] += coefficient * (new_y - old_y)

					q *= t
					r /= (1-t)

					old_x = new_x
					old_y = new_y

				BinomialCache = null

		Curvature(t)
			// This returns the unsigned scalar curvature on the curve at time t. This is computed
			// as
			//
			//   k(t) = ||Bézier'(t) x Bézier''(t)|| / ||Bézier'(t)||**3
			//
			// where x denotes the cross product. Even though this is a plane curve, for the
			// purposes of this method, when computing the cross product we treat the points
			// as though they are embedded on the Z-plane of an XYZ-space. That is, the X and
			// Y components can vary as needed but Z is fixed to zero.

			var
				// First derivative of Bézier(x) with respect to x, at time t.
				dx = 0
				dy = 0

				// Second derivative of Bézier(x) with respect to x, at time t.
				dx2 = 0
				dy2 = 0

				list/Data = TransformPointCache ? TransformPointCache : Points
				N = Length()

			if(N == 2)
				// If N == 2 (N < 2 is not allowed) then the Bézier curve is simply a line,
				// and thus has curvature zero.
				return 0

			if(t == 0)
				// Simplifications we can make when t == 0.

				dx = (N-1) * (Data[3] - Data[1])
				dy = (N-1) * (Data[4] - Data[2])

				// Bézier''(x) can be simplified to this.
				dx2 = (N-1)*(N-2) * (Data[5] - 2*Data[3] + Data[1])
				dy2 = (N-1)*(N-2) * (Data[6] - 2*Data[4] + Data[2])

			else if(t == 1)
				// Similar to above, but when t == 1.

				dx = (N-1) * (Data[Data.len-1] - Data[Data.len-3])
				dy = (N-1) * (Data[Data.len  ] - Data[Data.len-2])

				dx2 = (N-1)*(N-2) * (Data[Data.len-1] - 2*Data[Data.len-3] + Data[Data.len-5])
				dy2 = (N-1)*(N-2) * (Data[Data.len  ] - 2*Data[Data.len-2] + Data[Data.len-4])

			else
				// For the rest, we'll compute the first and second derivatives using the following
				// formulas, where b(n,k) = Choose(n,k) t**k * (1-t)**(n-k) and P(n) denotes the nth
				// control point.
				//
				//   Bézier'(t)  = n * Sum( b(n-1,k) * (P(k+1) - P(k)) ) from k = 0 to n-1.
				//   Bézier''(t) = n*(n-1) * Sum( b(n-2,k) * (P(k+2) - 2*P(k+1) + P(k)) ) from k = 0 to n-2.

				BinomialCache = new
				BinomialCache.len = N

				var
					dt_coefficient
					dt2_coefficient

					// Coordinates for the point P(k).
					P0_x = Data[1]
					P0_y = Data[2]

					// Coordinates for the point P(k+1).
					P1_x = Data[3]
					P1_y = Data[4]

					// Coordinates for the point P(k+2)
					P2_x = Data[5]
					P2_y = Data[6]

					q = 1
					r = (1-t)**(N-2)

				for(var/k = 0, k < (N-1), k ++)
					// Computing the sum for Bézier'(t)

					dt_coefficient = (N-1) * _Choose(N-2,k) * q*r
					dx += dt_coefficient * (P1_x - P0_x)
					dy += dt_coefficient * (P1_y - P0_y)

					// This must be done before Bézier''(t), because it has the value that
					// Bezier'(t) will have on the next iteration.
					r /= (1-t)

					if(k < (N-2))
						// And for Bézier''(t).

						dt2_coefficient = (N-1)*(N-2) * _Choose(N-3,k) * q*r
						dx2 += dt2_coefficient * (P2_x - 2*P1_x + P0_x)
						dy2 += dt2_coefficient * (P2_y - 2*P1_y + P0_y)

					// Shift P0, P1, and P2 appropriately.

					P0_x = P1_x
					P0_y = P1_y

					P1_x = P2_x
					P1_y = P2_y

					if(k < (N-3))
						P2_x = Data[2*k+7]
						P2_y = Data[2*k+8]

					// Both Bézier' and Bézier'' share the same value for q at each step, so this
					// must be done last.
					q *= t

				BinomialCache = null

			// Because A x B is perpendicular to both A and B, it's guaranteed that Bézier' x Bézier''
			// will lie entirely on the Z axis, as Bézier(t) lies entirely on the XY-plane. Thus, we
			// can simply compute the Z coordinate and take its absolute value in order to get
			// its length.
			. = abs(dx*dy2 - dx2*dy)

			// And this is |Bézier'(t)|.
			var/denom = sqrt(dx*dx + dy*dy)
			return . / (denom*denom*denom)