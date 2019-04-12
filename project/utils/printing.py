import time


def StrTimePrint(string, etime):

    out = ("{:60}").format(string), " ==>  ",
        ("{:4f}").format(etime), ("{:9}").format("   time: "),
        ("{:4}").format(time.strftime('%X %x %Z')[0:5])
    print (out)


def PrintStr(string, indent, boarder, preLines="", postLines=""):

    AppendedString = ""
    for i in range(indent):
        AppendedString = " " + AppendedString

    StringPrint = AppendedString + string
    if indent < 2:
        Outerboarder = ""
        if indent == 1:
            for i in range(int(np.ceil(len(string) / len(boarder))) + indent + 1):
                Outerboarder += boarder
        else:
            for i in range(int(np.ceil(len(string) / len(boarder))) + indent):
                Outerboarder += boarder

    else:
        AppendedString = ""
        for i in range(indent - 2):
            AppendedString = " " + AppendedString
        Outerboarder = AppendedString
        for i in range(int(np.ceil(len(string) / len(boarder))) + 4):
            Outerboarder += boarder
    print preLines + Outerboarder
    print StringPrint
    print Outerboarder + postLines
