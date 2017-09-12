import unittest
from util import parse_lastcol, ParseError

class UtilTest(unittest.TestCase):
    def test_parse_lastcol(self):
        good_string = 'a=b;c=d;e=12345;f="alphabet"'
        iffy_string = 'a=b;c=d;e=12345;f="alphabet";'
        
        good_parsed = {'a':'b', 'c':'d', 'e':'12345', 'f': '"alphabet"'}
        self.assertEqual(parse_lastcol(good_string), good_parsed)
        self.assertEqual(parse_lastcol(iffy_string), good_parsed)

        bad_string = 'a=b,c="d";'

        with self.assertRaises(ParseError):
            parse_lastcol(bad_string)


if __name__ == '__main__':
    unittest.main()
