from collections import defaultdict
import numpy as np
class Solution:
    def __init__(self):
        self.test =  [np.random.randn()*10 for i in range(100)]
        print(self.test)
    def hIndex(self, citations: list[int]) -> int:
        citations.sort(reverse=True)
        print(citations)
        maxi=0
        prev=float("-inf")
        print(citations)
        l=0
        r=len(citations)-1
        while l<=r:
            m=(r+l)//2
            print(m)
            if(citations[m]>=m+1 and (m+1==len(citations) or not citations[m+1]>=m+2)):
                return m+1
            if citations[m]>=m+1:
                l=m+1
                maxi=l
            else:
                r=m-1
                maxi=r
        return 0

let = Solution()
let.hIndex(let.test)
di = {}

