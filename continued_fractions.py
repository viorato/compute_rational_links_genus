"""
Continued fractions.
"""

from decimal import Decimal
from fractions import Fraction


class CFraction(list):
    """
    A continued fraction, represented as a list of integer terms.
    """

    def __init__(self, value, maxterms=15, cutoff=1e-10):
        if isinstance(value, (int, float, Decimal)):
            value = Decimal(value)
            remainder = int(value)
            self.append(remainder)

            while len(self) < maxterms:
                value -= remainder
                if value > cutoff:
                    value = Decimal(1) / value
                    remainder = int(value)
                    self.append(remainder)
                else:
                    break
        elif isinstance(value, (list, tuple)):
            self.extend(value)
        else:
            raise ValueError("CFraction requires number or list")

    def fraction(self, terms=None):
        "Convert to a Fraction."

        if terms is None or terms >= len(self):
            terms = len(self) - 1

        
        while terms > 0 and self[terms] == 0:
            terms -= 2
            
        if terms == 0:
            return Fraction(self[0],1)
        
        frac = Fraction(1,self[terms])

        for t in reversed(self[1:terms]):
            frac = 1 / (frac + t)

        frac += self[0]
        return frac

    def __float__(self):
        return float(self.fraction())

    def __str__(self):
        return "[%s]" % ", ".join([str(x) for x in self])



