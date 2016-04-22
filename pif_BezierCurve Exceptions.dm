#if !defined(PIF_NOPREFIX_GENERAL) && !defined(PIF_NOPREFIX_BEZIERCURVE)
pif_BezierCurve
#else
BezierCurve
#endif

	GenericBezierException
		parent_type = /exception

		New(_file, _line)
			file = _file
			line = _line

	InvalidConstructorArgumentFormatException
#if !defined(PIF_NOPREFIX_GENERAL) && !defined(PIF_NOPREFIX_BEZIERCURVE)
		parent_type = /pif_BezierCurve/GenericBezierException
#else
		parent_type = /BezierCurve/GenericBezierException
#endif

		name = "Invalid Constructor Argument Format"
		desc = "Invalid format for constructor argument. Expected either a list or an argument list of coordinate pairs."

	InvalidArgumentFormatException
#if !defined(PIF_NOPREFIX_GENERAL) && !defined(PIF_NOPREFIX_BEZIERCURVE)
		parent_type = /pif_BezierCurve/GenericBezierException
#else
		parent_type = /BezierCurve/GenericBezierException
#endif

		name = "InvalidArgumentFormat"
		desc = "Invalid or unknown argument format for method."

	InvalidArgumentDataTooFewPointsException
#if !defined(PIF_NOPREFIX_GENERAL) && !defined(PIF_NOPREFIX_BEZIERCURVE)
		parent_type = /pif_BezierCurve/GenericBezierException
#else
		parent_type = /BezierCurve/GenericBezierException
#endif

		name = "Invalid Argument Data (Too Few Points)"
		desc = "Invalid argment data: Must provde at least two points (four coordinates)."

	InvalidArgumentDataOddNumberException
#if !defined(PIF_NOPREFIX_GENERAL) && !defined(PIF_NOPREFIX_BEZIERCURVE)
		parent_type = /pif_BezierCurve/GenericBezierException
#else
		parent_type = /BezierCurve/GenericBezierException
#endif

		name = "Invalid Argument Data (Odd Number Of Coordinates)"
		desc = "Invalid argment data: N points require 2N coordinates."