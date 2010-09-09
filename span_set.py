import abc
import collections
import itertools
import operator
import unittest

# A Span represents a contiguous range of text in a document, by
# storing its zero-based start and end points.  The end point is
# non-inclusive, so for example, if the document is the string
# "absurd", the span (2, 4) represents the substring "su".
class Span(collections.namedtuple('Span', 's e')):
    def __new__(cls, s, e):
        assert 0 <= s <= e
        return cls.__bases__[0].__new__(cls, s, e)


# A SpanRect represents the set of spans (s, e) such that smin <= s <=
# smax and emin <= e <= emax.  It is called a SpanRect because when
# represented in "span space" (the 2-dimensional space defined by the
# s and e axes), it takes the shape of a rectangle intersected with
# the region s <= e.
#
# A SpanRect is required to represent a non-empty set of spans.  This
# means that smin <= smax, emin <= emax, and smin <= emax.  In
# addition, to ensure uniqueness we require that 0 <= smin, smax <=
# emax, smin <= emin.
class SpanRect(collections.namedtuple('SpanRect', 'smin emin smax emax')):
    def __new__(cls, smin, emin, smax, emax):
        assert 0 <= smin <= smax <= emax
        assert smin <= emin <= emax
        return cls.__bases__[0].__new__(cls, smin, emin, smax, emax)

    # Figure out if a SpanRect contains the given Span.
    def __contains__(self, span):
        return self.smin <= span.s <= self.smax and self.emin <= span.e <= self.emax

    # Figure out if two SpanRects intersect.
    @staticmethod
    def intersects(r1, r2):
        smin = max(r1.smin, r2.smin)
        smax = min(r1.smax, r2.smax)
        if smin > smax:
            return False
        emin = max(r1.emin, r2.emin)
        emax = min(r1.emax, r2.emax)
        if emin > emax:
            return False
        return True

    # Return a list of SpanRects representing the set subtraction of
    # r1 and r2.
    @staticmethod
    def subtract(r1, r2):
        result = []
        # Figure out critical indices in the s and e dimensions where
        # the truth value (s, e) in (r1 - r2) might change.  This
        # divides span space into a set of nonoverlapping rectangles
        # such that each rectangle is either entirely included in the
        # result or entirely absent from the result.
        critical_s_points = sorted(set([r1.smin, r1.smax + 1, r2.smin, r2.smax + 1]))
        critical_e_points = sorted(set([r1.emin, r1.emax + 1, r2.emin, r2.emax + 1]))
        # Now test each of these rectangles at its upper left corner
        # to see whether the result should contain it.
        for i in xrange(len(critical_s_points) - 1):
            for j in xrange(len(critical_e_points) - 1):
                smin = critical_s_points[i]
                emax = critical_e_points[j + 1] - 1
                if smin <= emax:
                    span = Span(smin, emax)
                    if span in r1 and span not in r2:
                        smax = min(critical_s_points[i + 1] - 1, emax)
                        emin = max(critical_e_points[j], smin)
                        result.append(SpanRect(smin, emin, smax, emax))
        return result

# Abstract base class defining the interface for any representation of
# a set of Spans.
class SpanSet(object):
    __metaclass__ = abc.ABCMeta

    # Yield a stream of SpanRects whose union is self.
    #
    # Optional argument restrict_rect, informally, is a SpanRect
    # representing a subset of the SpanSet to be queried.  Formally,
    # query(self, restrict_rect) yields a stream of SpanRects whose
    # union, when intersected with restrict_rect, is equal to the
    # intersection of self and restrict_rect.
    @abc.abstractmethod
    def query(self, restrict_rect = None):
        pass


# Abstract base class for any SpanSet that can be mutated in place by
# adding and removing SpanRects.
class MutableSpanSet(SpanSet):

    # Insert a new SpanRect into the set, so that the set of spans
    # represented by the SpanRect is unioned with the set of spans
    # represented by self.
    @abc.abstractmethod
    def insert(self, r):
        pass

    # Delete a SpanRect from the set, so that the set of spans
    # represented by the SpanRect is subtracted from the set of spans
    # represented by self.
    @abc.abstractmethod
    def delete(self, r):
        pass

    # Make a copy of this set that can be modified independently of
    # self.
    @abc.abstractmethod
    def copy(self):
        pass


# An inefficient implementation of MutableSpanSet based on storing
# SpanRects in a simple list.  Probably too inefficient for practical
# use, but good for experimenting with.
class SpanRectList(MutableSpanSet):

    # SpanRectList() with no arguments creates an object representing
    # the empty set.  SpanRectList with a list (or other enumerable)
    # of SpanRects creates an object representing the union of those
    # rects.  This function makes a copy of the contents of its
    # argument.
    def __init__(self, rects = []):
        self.__rects = list(rects)

    def query(self, restrict_rect = None):
        for r in self.__rects:
            if restrict_rect and not SpanRect.intersects(restrict_rect, r):
                continue
            yield r

    def insert(self, r):
        self.__rects.append(r)

    def delete(self, r):
        rects = []
        for r2 in self.__rects:
            rects.extend(SpanRect.subtract(r2, r))
        self.__rects = rects

    def copy(self):
        return SpanRectList(self.__rects)


class SpanRectTestHelper(object):
    @staticmethod
    def all_spans(size):
        result = []
        for s in xrange(size):
            for e in xrange(s, size):
                result.append(Span(s, e))
        return result

    @staticmethod
    def all_rects(size):
        result = []
        for smin in xrange(size):
            for emin in xrange(smin, size):
                for smax in xrange(smin, size):
                    for emax in xrange(max(smax, emin), size):
                        result.append(SpanRect(smin, emin, smax, emax))
        return result


class SpanRectTests(SpanRectTestHelper, unittest.TestCase):
    def test_contains(self):
        r = SpanRect(2, 4, 6, 8)
        for s, e in ((2, 4), (4, 4), (4, 6), (6, 6), (2, 8), (6, 8)):
            self.assertTrue(Span(s, e) in r, repr((s, e)))
        for s, e in ((3, 3), (1, 6), (7, 7), (4, 9)):
            self.assertFalse(Span(s, e) in r, repr((s, e)))

    @classmethod
    def slow_intersects(cls, r1, r2, size):
        for span in cls.all_spans(size):
            if span in r1 and span in r2:
                return True
        return False
        
    def test_intersects(self):
        size = 5
        rs = self.all_rects(size)
        for i in xrange(len(rs)):
            for j in xrange(i, len(rs)):
                r1 = rs[i]
                r2 = rs[j]
                expected = self.slow_intersects(r1, r2, size)
                actual = SpanRect.intersects(r1, r2)
                self.assertEqual(
                    expected, actual,
                    'SpanRect.intersects({0!r}, {1!r}) == {2}, expected {3}'.format(r1, r2, actual, expected))
                actual = SpanRect.intersects(r2, r1)
                self.assertEqual(
                    expected, actual,
                    'SpanRect.intersects({0!r}, {1!r}) == {2}, expected {3}'.format(r2, r1, actual, expected))

    def test_subtract(self):
        size = 5
        rs = self.all_rects(size)
        ss = self.all_spans(size)
        for i in xrange(len(rs)):
            for j in xrange(i, len(rs)):
                r1 = rs[i]
                r2 = rs[j]
                result = SpanRect.subtract(r1, r2)
                for span in ss:
                    self.assertEqual(span in r1 and span not in r2, any(span in r for r in result))


# Base class for testing any MutableSpanSet-derived object.
class MutableSpanSetTests(SpanRectTestHelper):
    __metaclass__ = abc.ABCMeta

    # Override in derived class to call the appropriate constructor
    # and return an empty set.
    @abc.abstractmethod
    def mk_set(self):
        pass

    def check_query_result(self, size, ss, f, context, restrict = None):
        # Run a query against ss, and verify that the union of the
        # resulting rects is equivalent to the set defined by the
        # membership test f.  Check all spans up to the given size.
        query_result = list(ss.query(restrict))
        if restrict is not None:
            for r in query_result:
                self.assertTrue(SpanRect.intersects(r, restrict))
        for span in self.all_spans(size):
            if restrict is None or span in restrict:
                self.assertEqual(
                    any(span in r for r in query_result), f(span),
                    ("{0} produces {1!r}, which, when queried with restriction {2!r}, produces {3!r}.  "
                     "Mismatch at {4!r}.").format(
                        context, ss, restrict, query_result, span))

    def insert_delete_n_test(self, size, n):
        for rs in itertools.product(self.all_rects(size), repeat=n):
            ss = self.mk_set()
            for i, r in enumerate(rs):
                if i % 2 == 0:
                    ss.insert(r)
                else:
                    ss.delete(r)
            def f(span):
                in_set = False
                for i, r in enumerate(rs):
                    if span in r:
                        in_set = (i % 2 == 0)
                return in_set
            context = 'insert/delete {0!r}'.format(rs)
            self.check_query_result(size, ss, f, context)
            for restrict in rs:
                self.check_query_result(size, ss, f, context, restrict)

    def test_insert_delete_1(self):
        self.insert_delete_n_test(8, 1)

    def test_insert_delete_2(self):
        self.insert_delete_n_test(3, 2)

    def test_insert_delete_3(self):
        self.insert_delete_n_test(2, 3)

    def test_insert_delete_4(self):
        self.insert_delete_n_test(2, 4)


class SpanRectListTests(MutableSpanSetTests, unittest.TestCase):
    def mk_set(self):
        return SpanRectList()


if __name__ == '__main__':
    unittest.main()
