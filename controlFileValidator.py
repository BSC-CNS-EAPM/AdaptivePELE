import validatorBlockNames
import blockNames
import json

def validate(control_file):
    
    with open(control_file,'r') as f:
        jsonFile = f.read()
    parsedJSON = json.loads(jsonFile)
    for block in dir(validatorBlockNames.ControlFileParams):
        if block.startswith('__'):
            continue
        block_obj = getattr(validatorBlockNames,
                            eval("validatorBlockNames.ControlFileParams.%s"%block))
        controlfile_obj = eval("blockNames.ControlFileParams.%s"%block)
        #try:
        if block == "generalParams":
            validateGeneralBlock(block_obj, parsedJSON[controlfile_obj])
        else:
            validateBlock(block_obj, parsedJSON[controlfile_obj])
        #except KeyError:
        #    raise KeyError("Block %s not found in control file!"%controlfile_obj)
    print "Congratulations! No errors found!"
    return True

def validateBlock(blockName, controlFileBlock):
    blockType = controlFileBlock["type"]
    #Check if type selected is valid
    try:
        assert isinstance(blockType,unicode),"Type for %s should be %s and instead is %s"%(blockType, 'unicode', type(blockType).__name__)
    except KeyError:
        raise KeyError("Type %s in %s not found."%(blockType,blockName.__name__))
    # check for mandatory parameters
    for mandatory,value in blockName.types[blockType].iteritems():
        try:
            assert isinstance(controlFileBlock['params'][mandatory],eval(value)),"Type for %s should be %s and instead is %s"%(mandatory, value, type(controlFileBlock['params'][mandatory]).__name__)
        except KeyError as err:
            raise KeyError("%s missing: Mandatory parameter %s in %s not found."%(err.message,mandatory,blockName.__name__)) 
    # check rest of parameters specified
    for param, value in controlFileBlock["params"].iteritems():
        try:
            assert isinstance(value,eval(blockName.params[param])),"Type for %s should be %s and instead is %s"%(param, blockName.params[param], type(value).__name__)
        except KeyError:
            raise KeyError("Parameter %s not recognized."%param) 

    for block in dir(blockName):
        if not block.startswith('__') and block != "params" and block != "types":
            types_dict = eval("blockName.%s"%block)["types"]
            params_dict = eval("blockName.%s"%block)["params"]
            try:
                blockType = controlFileBlock[block]["type"]
                try:
                    assert isinstance(blockType,str),"Type for %s should be %s and instead is %s"%(blockType, 'str', type(blockType).__name__)
                    typecheck = types_dict[blockType]
                except KeyError:
                    raise KeyError("Type %s in %s not found."%(blockType,blockName.__name__))
            # check rest of parameters specified
                for param, value in controlFileBlock[block].iteritems():
                    try:
                        assert isinstance(value,eval(params_dict[param])),"Type for %s should be %s and instead is %s"%(param, params_dict[param], type(value).__name__)
                    except KeyError:
                        raise KeyError("Parameter %s not recognized."%param) 
            except KeyError:
                #The parameters blocks for density and threshold calculator are
                #not mandatory
                pass

def validateGeneralBlock(blockName, controlFileBlock):
    for key,value in controlFileBlock.iteritems():
        try:
            assert isinstance(controlFileBlock[key],eval(blockName.params[key])),"Type for %s should be %s and instead is %s"%(key, value, type(controlFileBlock[key]).__name__)
        except KeyError:
            raise KeyError("Mandatory parameter %s in GeneralParams not found."%key) 

if __name__ == "__main__":
    controlFiles = ["tests/data/3ptb_data/integrationTest%i.conf"%i for i in range(1,4)]
    for contfile in controlFiles:
        print "Validating control file %s"%contfile
        validate(contfile)
