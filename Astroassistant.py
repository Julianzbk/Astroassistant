"""
    Astroassistant2 in Python, v1.0
    A module used in calculations for Astrophysics, Rocketry, and Orbital Mechanics.
    Not for precise scientific usage

    Written by Julianzbk, github.com/Julianzbk

    Call functions for individual usage, run main() for a complete user interface.
"""
import math
import csv

class KwargsError(Exception):
    """An invalid keyword argument was entered into a function.
    """

def get_constants(*args:str) -> tuple[int,int] or None:
    """Internal function, reads from constant.csv and returns constants.
    Args (Optional):
        *args: Name of the primary body to search for.
    If no Args:
        Prompts the user for the name of the primary.
    Returns:
        constants: A tuple of the body's Standard Gravitaional Parameter (μ = G*M) and radius (m).
        None: If args was entered but is empty.
    Raises:
        NameError: The name provided was not found in file.
    """
    constants = ()
    while not constants:
        if not args:
            name = input('Enter name of the primary body:\n')
        else:
            name = args[0]
            if name == '':
                return None

        try:
            with open('constants.csv', 'r', encoding="utf-8") as file:
                lines = csv.reader(file, quotechar='"', quoting=csv.QUOTE_ALL, skipinitialspace=True, delimiter=',')
                for i, row in enumerate(lines):
                    for value in row[:-1]:
                        if value == name:
                            constants = tuple([int(float(num)) for num in row[-1].split(',')]) # Store (G*M, Radius)
                            if len(constants) != 2:
                                raise ValueError('Invalid constant tuple.')

        except (TypeError, SyntaxError, ValueError):
            print('ERROR: Check constants.csv for formatting issues.')
            print(f'Traceback (Line {i + 1}):\n {row}')
            input('Press Enter to proceed:')
            raise
        except FileNotFoundError as fnf:
            raise FileNotFoundError('constants.csv not found, please download file and place it under the same directory.') from fnf

        if not constants:
            if args:
                raise NameError('Body name not found, try again.')
            print('Body name not found, try again.')
    return constants

def kwargs_handle(kwargs:dict, accepted_kwords:list) -> dict:
    """Internal Function. Checks provided keyword arguments against a list of accepted keywords, then returns a formatted dict of kwargs.\n
    Args:
        kwargs: A dictionary of keyword arguments.
        accepted_kwords: A list of accepted keywords for a function.
    Returns:
        new_kwargs: A dictionary of formatted keyword arguments.
    Raises:
        KwargsError: The keywords do not match exactly to accepted_kwords.
    """
    accepted_kwords = sorted(list({'ui'}.union(set(accepted_kwords))))

    kwords = []
    new_kwargs = {}
    for k in kwargs:
        value = kwargs[k]
        k = k.lower()
        kwords.append(k)
        new_kwargs.update({k: value})
    if 'ui' not in new_kwargs:
        new_kwargs.update({'ui':False})
        kwords.append('ui')

    if sorted(kwords) != accepted_kwords:
        raise KwargsError('Invalid keyword arguments entered.')
    return new_kwargs

def conjunct_check(kwords:list, kwargs:dict) -> bool or int:
    """Internal Function. Checks if the keywords provided exists in conjunction within the provided kwargs dictionary.
    Args:
        kwords: A list of keywords that must be used in conjuction.
        kwargs: The dict of keyword arguments passed to the original function.
    Returns:
        True: If all kwords are found within kwargs.
        False: If some but not all kwords are found within kwargs.
        -1: If no kwords are found within kwargs.
    """
    ANDcond = True
    ORcond = False
    for k in kwords:
        ANDcond = ANDcond and (k in kwargs)
        ORcond = ORcond or (k in kwargs)

    if ANDcond:
        return True
    elif ORcond:
        return False
    else:
        return -1

def unit_format(quantity:str) -> float:
    """Internal Function. Expands quantities of distance with a unit.\n
    Args:
        quantity: The quantity to convert, with a unit as its last 2 characters.
    Returns:
        The expanded quantity in meters.
    """
    if not isinstance(quantity, str):
        return quantity
    
    if quantity[0] == "'" and quantity[-1] == "'":
        quantity = quantity[1:-1]

    if quantity[-2:] == 'km':
        quantity = quantity[0:-2]
        if quantity.isdigit():
            return int(quantity)*1000
        return float(quantity)*1000

    if quantity[-2:] == 'Mm' or quantity[-2:] == 'mm':
        quantity = quantity[0:-2]
        if quantity.isdigit():
            return int(quantity)*1000000
        return float(quantity)*1000000

    if quantity[-2:] == 'AU' or quantity[-2:] == 'au':
        quantity = quantity[0:-2]
        if quantity.isdigit():
            return int(quantity)*149_597_870_700
        return float(quantity)*149_597_870_700

    if quantity.isdigit():
        return int(quantity)
    else:
        return float(quantity)

def time_format(timeSec:float) -> str:
    """Internal Function. Formats durations of time into a readable string.\n
    Args:
        timeSec: The duration of time in seconds.
    Returns:
        A string that represents a duration of time in years/days/hours/minutes/seconds(1 decimal places).
    """
    if timeSec <= 3600:
        timeMin = timeSec // 60
        timeSec = timeSec % 60
        return f'{int(timeMin)} min {timeSec:.1f} sec'
    elif timeSec <= 86400:
        timeHour = timeSec // 3600
        timeMin = timeSec % 3600 // 60
        timeSec = timeSec % 3600 % 60
        return f'{int(timeHour)} hr {int(timeMin)} min {timeSec:.1f} sec'
    elif timeSec <= 31536000:
        timeDay = timeSec // 86400
        timeHour = timeSec % 86400 // 3600
        timeMin = timeSec % 86400 % 3600 // 60
        timeSec = timeSec % 86400 % 3600 % 60
        return f'{int(timeDay)} day {int(timeHour)} hr {int(timeMin)} min {timeSec:.1f} sec'
    else:
        timeYear = timeSec // 31536000
        timeDay = timeSec % 31536000 // 86400
        timeHour = timeSec % 31536000 % 86400 // 3600
        timeMin = timeSec % 31536000 % 86400 % 3600 // 60
        timeSec = timeSec % 31536000 % 86400 % 3600 % 60
        return f'{int(timeYear)} year {int(timeDay)} day {int(timeHour)} hr {int(timeMin)} min {timeSec:.1f} sec'


def orbital_velocity(**kwargs) -> float:
    """Calculates the current orbital velocity of a vessel.\n
    Args (Keyword, Optional):
        r: The current altitude of the vessel (Above Sea Level), aka
           the radius of its orbit at this point (m).
        a: The Semi-Major Axis of the vessel's orbit (Above Sea Level), aka
           the average altitude of the vessel (m).
        body: The name of the astronomical body that the vessel orbits.
    Can be used individually:
        center (optional): Whether or not 'r' and 'a' are entered as altitudes Above Sea Level. (default False)
    For orbit around a custom primary body (Cannot be used in conjunction with 'body'):
        mu: The Standard Gravitational Parameter (μ = G*M) of the primary body. Must be used in conjunction with:
        bodyr: The radius of the primary body (m). (Obsolete if center=True)\n
    Returns:
        v: The current Orbital Velocity of the vessel (m/s).
    """
    if kwargs: # You should probably just fold this up...
        # Check for optional keywords.
        conjunct_kwords = ['mu','bodyr']
        conjunct = conjunct_check(conjunct_kwords, kwargs)
        if conjunct is True:
            constants = (unit_format(kwargs['mu']), unit_format(kwargs['bodyr']))
            del kwargs['mu'], kwargs['bodyr']
            accepted_kwords = ['r','a']
        elif conjunct is False:
            raise KwargsError('Invalid keyword arguments entered. "mu" and "bodyr" must be used in conjunction.')
        else:
            accepted_kwords = ['r','a','body']
            constants = ()

        if 'center' in kwargs:
            center = kwargs['center']
            if not isinstance(center, bool):
                raise KwargsError('Invalid keyword arguments entered: "center" is not bool')
            del kwargs['center']
        else:
            center = False
            if not accepted_kwords:
                accepted_kwords = ['r','a','body']

        # Check for center only.
        conjunct = conjunct_check(accepted_kwords, kwargs)
        if conjunct == -1:
            kwargs = False # !!
            ui = True
            constants = get_constants()
        else: # Check for required keywords.
            kwargs = kwargs_handle(kwargs, accepted_kwords)
            r = kwargs['r']
            a = kwargs['a']
            ui = kwargs['ui']
            if not constants:
                constants = get_constants(kwargs['body'])
    else:
        center = False
        ui = True
        constants = get_constants()

    if not kwargs:
        if center:
            r = input('Enter current distance from center: ')
            a = input('Enter Semi-Major Axis from center:  ')
        else:
            r = input('Enter current altitude ASL: ')
            a = input('Enter Semi-Major Axis ASL:  ')
    r = unit_format(r)
    a = unit_format(a)
    if not center:
        r += constants[1]
        a += constants[1]

    v = math.sqrt(constants[0]*(2/r - 1/a))
    if ui:
        print(f'\nOrbital Velocity: {v:.2f} m/s')
    return v

def inclination_dv(**kwargs) -> float:
    """Calculates the Delta-V required to perform a change of inclination maneuver for a vessel.\n
    Args (Keyword, Optional):
        vi: The initial orbital velocity of the vessel (m/s).
        vf: The final orbital velocity of the vessel (m/s). (None = vi)
        inc: The change in inclination of the maneuver (°).
    Returns:
        v: The Delta-V required to perform the inclination maneuver. (m/s)
    Optional Kwargs (Cannot be used in conjunction with 'vi' and 'vf'):
        ri: The current altitude of the vessel (Above Sea Level, m).
        ai: The Semi-Major Axis of the current orbit (Above Sea Level, m).
        rf (Optional): The final altitude of the vessel (Above Sea Level, m). (default ri)
        af (Optional): The Semi-Major Axis of the final orbit (Above Sea Level, m). (default ai)
        body: The name of the astronomical body that the vessel orbits.
    """
    if kwargs:
        conjunct_kwords = ['ri','ai','body']
        conjunct = conjunct_check(conjunct_kwords, kwargs)
        if conjunct:
            vi = orbital_velocity(r=kwargs['ri'], a=kwargs['ai'], body=kwargs['body'])

            conjunct = conjunct_check(['rf','af'], kwargs)
            if conjunct:
                vf = orbital_velocity(r=kwargs['rf'], a=kwargs['af'], body=kwargs['body'])
            elif not conjunct:
                raise KwargsError('Invalid keyword arguments entered. "rf" and "af" must be used in conjunction.')
            else:
                vf = vi

            if 'ui' not in kwargs: # Seperate logic for default ui value.
                kwargs['ui'] = False

        elif not conjunct:
            raise KwargsError('Invalid keyword arguments entered. "r-", "a-", and "body" must be used in conjunction.')
        else:
            vi = False
            accepted_kwords = ['vi','vf','inc']
            kwargs = kwargs_handle(kwargs, accepted_kwords)
            if not vi:
                vi = kwargs['vi']
                if kwargs['vf'] is None:
                    vf = vi
                else:
                    vf = kwargs['vf']

        inc = kwargs['inc']
        ui = kwargs['ui']
    else:
        ui = True

    if not kwargs:
        vi = float(input('Enter velocity of spacecraft before maneuver (m/s):\n'))
        vf = float(input('Enter velocity of spacecraft after maneuver (m/s):\n'))
        inc = float(input('Enter change of inclination in degrees:\n'))
    v = math.sqrt(vi**2 + vf**2 - 2*vi*vf*math.cos(math.radians(inc)))
    if ui:
        print(f'\nDelta-V of maneuver: {v:.6g} m/s')
    return v

def orbital_period(**kwargs) -> float:
    """Calculates the Orbital Period of a vessel.
    Args (Keyword, Optional):
        a: The Semi-Major Axis of the vessel's orbit (Above Sea Level), aka
           the average altitude of the vessel (m).
        body: The name of the astronomical body that the vessel orbits.
    Returns:
        timeSec: Orbital Period in seconds, call time_format() to return a formatted string.
    Can be used individually:
        center (optional): Whether or not 'a' is entered as altitude Above Sea Level. (default False)\n
    For orbit around a custom primary body (Cannot be used in conjunction with 'body'):
        mu: The Standard Gravitational Parameter (μ = G*M) of the primary body. Must be used in conjunction with:
        bodyr: The radius of the primary body (m). (Obsolete if center=True)
    """
    if kwargs: # You should probably just fold this up...
        # Check for optional keywords.
        conjunct_kwords = ['mu','bodyr']
        conjunct = conjunct_check(conjunct_kwords, kwargs)
        if conjunct is True:
            constants = (unit_format(kwargs['mu']), unit_format(kwargs['bodyr']))
            del kwargs['mu'], kwargs['bodyr']
            accepted_kwords = ['a']
        elif conjunct is False:
            raise KwargsError('Invalid keyword arguments entered. "mu" and "bodyr" must be used in conjunction.')
        else:
            accepted_kwords = ['a','body']
            constants = ()

        if 'center' in kwargs:
            center = kwargs['center']
            if not isinstance(center, bool):
                raise KwargsError('Invalid keyword arguments entered: "center" is not bool')
            del kwargs['center']
        else:
            center = False
            if not accepted_kwords:
                accepted_kwords = ['a','body']

        # Check for center only.
        conjunct = conjunct_check(accepted_kwords, kwargs)
        if conjunct == -1:
            kwargs = False
            ui = True
            constants = get_constants()
        else: # Check for required keywords.
            kwargs = kwargs_handle(kwargs, accepted_kwords)
            a = kwargs['a']
            ui = kwargs['ui']
            if not constants:
                constants = get_constants(kwargs['body'])
    else:
        center = False
        ui = True
        constants = get_constants()

    if not kwargs:
        if center:
            a = input('Enter Semi-Major Axis from center: ')
        else:
            a = input('Enter Semi-Major Axis ASL: ')
    a = unit_format(a)
    if not center:
        a += constants[1]

    timeSec = math.tau * math.sqrt(a**3/constants[0])
    if ui:
        print('\nOrbital Period:', time_format(timeSec))
    return timeSec

def two_body(**kwargs) -> float:
    """Calculates the gravitational Force between two bodies of mass.
    Args (Keyword, Optional):
        m1: The Mass (kg) or Name of the first object.
        m2: The Mass (kg) or Name of the second object.
        r: The distance between the objects (center-center) (m).
    Returns:
        f: The gravitational Force exerted by the objects on each other (N).
    """
    if kwargs:
        accepted_kwords = ['m1','m2','r']
        kwargs = kwargs_handle(kwargs, accepted_kwords)
        m1 = str(kwargs['m1'])
        m2 = str(kwargs['m2'])
        r = unit_format(kwargs['r'])
        ui = kwargs['ui']
    else:
        ui = True

    if not kwargs:
        m1 = input('Enter Mass (kg) or Name of first object.\n')
        m2 = input('Enter Mass (kg) or Name of second object.\n')
        r = input('Enter distance between objects.\n')

    if m1.replace('.','').isnumeric():
        m1 = float(m1)
    else:
        m1 = get_constants(m1)[0] / 6.6743e-11

    if m2.replace('.','').isnumeric():
        m2 = float(m2)
    else:
        m2 = get_constants(m2)[0] / 6.6743e-11
    r = unit_format(r)

    f = (6.6743e-11*m1*m2)/(r**2)
    if ui:
        print(f'\nGravitational Force: {f:.6g} N\n')
    return f

def rocket_equation(**kwargs) -> float:
    """Calculates the Delta-V available to a rocket stage. 
    Args (Keyword, Optional):
        wetmass: The Mass of the vehicle when fully loaded with fuel (kg).
        drymass: The Mass of the vehicle when emptied of fuel (kg).
        isp: The Specific Impulse of the rocket engine (N*s/kg).
            OR
        ve: The Exhaust Velocity of the rocket engine (m/s).
    Returns:
        v: The stage Delta-V (m/s).
    """
    if kwargs:
        conjunct_kwords = ['isp','ve']
        conjunct = conjunct_check(conjunct_kwords, kwargs)
        if not conjunct:
            if 'isp' in kwargs:
                isp = kwargs['isp']
                para = 'isp'
                del kwargs['isp']
            if 've' in kwargs:
                ve = kwargs['ve']
                para = 've'
                del kwargs['ve']
        elif conjunct:
            raise KwargsError('Invalid keyword arguments entered. "isp" and "ve" cannot both be used.')
        else:
            raise KwargsError('Invalid keyword arguments entered. Must have either "isp" or "ve".')
        
        accepted_kwords = ['wetmass','drymass']
        kwargs = kwargs_handle(kwargs, accepted_kwords)
        wetmass = kwargs['wetmass']
        drymass = kwargs['drymass']
        ui = kwargs['ui']
    else:
        ui = True

    if not kwargs:
        wetmass = float(input('Enter vehicle wet mass (kg): '))
        drymass = float(input('Enter vehicle dry mass (kg): '))
        para = input('Enter Specific Impulse or Exhaust Velocity, in the form isp=... ve=...\n')

    if para[:3] == 'isp':
        if not kwargs:
            isp = float(para.split('=')[-1])
        v = isp * 9.80665 * math.log(wetmass/drymass)
    elif para[:2] == 've':
        if not kwargs:
            ve = float(para.split('=')[-1])
        v = ve * math.log(wetmass/drymass)
    else:
        raise ValueError('Improper input.')

    if ui:
        print(f'\nStage Delta-V: {v:.2f} m/s')
    return v

def orbital_eccentricity(**kwargs) -> float:
    """Calculates the Orbital Eccentricity of a vessel.
    Args (Keyword, Optional):
        ap: The Apoapsis (farthest approach) of the orbit. (ASL/center)
        pe: The Periapsis (closest approach) of the orbit. (ASL/center)
        body (Optional): The name of the astronomical body that the vessel orbits. (center if empty)
    Returns:
        ecc: The Eccentricity of the orbit.
    """
    if kwargs:
        if 'body' in kwargs:
            constants = get_constants(kwargs['body'])
            del kwargs['body']
        else:
            constants = None
        accepted_kwords = ['ap','pe']
        kwargs = kwargs_handle(kwargs, accepted_kwords)
        ap = kwargs['ap']
        pe = kwargs['pe']
        ui = kwargs['ui']
    else:
        ui = True
        constants = get_constants(input('Enter Body name, leave blank if irrelevant.\n'))

    if constants is None and not kwargs:
        ap = input('Enter Apoapsis from center:  ')
        pe = input('Enter Periapsis from center: ')
    elif not kwargs:
        ap = input('Enter Apoapsis ASL:  ')
        pe = input('Enter Periapsis ASL: ')

    ap = unit_format(ap)
    pe = unit_format(pe)
    if constants:
        ap += constants[1]
        pe += constants[1]

    ecc = (ap - pe)/(ap + pe)
    if ui:
        print('\nEccentricity: ', end='')
        print(f'{ecc:.6f}')
    return ecc

def hyperbolic_velocity(**kwargs) -> float:
    """Calculates the Hyperbolic Excess Velocity of a vessel.
    Args (Keyword, Optional):
        v: The current orbital velocity of the vessel (m/s).
        r: The current altitude of the vessel (Above Sea Level), aka
           the radius of its orbit at this point (m). 
        body: The name of the astronomical body that the vessel orbits.
        center (Optional): Whether or not 'r' is entered as altitude Above Sea Level. (default False)
    Returns:
        vinf: The Hyperbolic Excess Velocity of the vessel.
        -1.0: If the vessel is below escape velocity.
    """
    if kwargs:
        if 'center' in kwargs:
            center = kwargs['center']
            if not isinstance(center, bool):
                raise KwargsError('Invalid keyword arguments entered: "center" is not bool')
            del kwargs['center']
        else:
            center = False

        accepted_kwords = ['v','r','body']
        kwargs = kwargs_handle(kwargs, accepted_kwords)
        v = kwargs['v']
        r = kwargs['r']
        ui = kwargs['ui']
        constants = get_constants(kwargs['body'])
    else:
        constants = get_constants()
        center = False
        ui = True
    
    if not kwargs:
        v = float(input('Enter current Velocity: '))
        r = input('Enter current altitude ASL: ')
    r = unit_format(r)
    if not center:
        r += constants[-1]

    vesc = orbital_velocity(r=r, a=r, mu=constants[0], bodyr=constants[1], center=True) * math.sqrt(2)
    if v < vesc:
        if ui:
            print('\nCurrent Velocity is less than Escape Velocity.')
        return -1.0
    vinf = math.sqrt(v**2 - vesc**2)
    if ui:
        print(f'\nHyperbolic Excess Velocity: {vinf:.6g} m/s')
    return vinf

def create_body(**kwargs) -> str:
    '''Creates a custom Primary Body and stores it in constants.csv.
    Args (Keyword, Optional):
        name: The name of the body.
        r: The radius of the body. (m)
        mass: The mass of the body. (kg)
            OR
        mu: The Standard Gravitional Parameter of the body. (kg^3*m*s^-2)
    Returns:
        line: The string that was written into constants.csv.
    '''
    if kwargs:
        # Disjunct check for optional kwords.
        conjunct = conjunct_check(['mass','mu'], kwargs)
        if not conjunct:
            if 'mass' in kwargs:
                m = f"{kwargs['mass'] * 6.67423e-11:e}"
                del kwargs['mass']
            elif 'mu' in kwargs:
                m = f"{kwargs['mu']:e}"
                del kwargs['mu']
        else:
            raise KwargsError('Invalid keyword arguments entered. "mass" and "mu" cannot both be used.')
        accepted_kwords = ['name','r']
        kwargs_handle(kwargs, accepted_kwords)
        names = kwargs['name']
        r = kwargs['r']
        ui = kwargs['ui']
    else:
        ui = True

    if not kwargs:
        names = input('Enter Name of the body, or Names seperated by a comma:\n')
        m = input('\nEnter Mass or Standard Gravitational Parameter of the body, in the form: mass=... mu=...\n')
        if m[:4] == 'mass':
            m = f"{float(m.split('=')[-1]) * 6.67423e-11:.7e}"
        elif m[:2] == 'mu':
            m = f"{float(m.split('=')[-1]):.7e}"
        else:
            raise ValueError('Improper input.')
        r = input('\nEnter Radius of the body.\n')

    names = names.split(',')
    line = []
    for name in names:
        line.append(name.strip() + ', ')
    line += ['"', m, ', ', f"{unit_format(r):.3e}", '"']
    line = ''.join(line)

    try:
        with open('constants.csv', 'a', encoding="utf-8") as file:
            from os import getlogin
            from time import ctime
            file.write('\n\n')
            file.write(f'# Custom Body created by {getlogin()} on ({ctime()}):\n')
            file.write(line)
            file.write('\n')
    except FileNotFoundError as fnf:
        raise FileNotFoundError('constants.csv not found, please download file and place it under the same directory.') from fnf

    if ui:
        print('\nCustom Body created:')
        print(line)
    return line


def main():
    """
        Main function for running Astroassistant with an UI.
    """
    print('Welcome to Astroassistant, an astrophysics calculator.\n')
    functions = [(1, 'Orbital Velocity', orbital_velocity), 
                 (2, 'Inclination Delta-V', inclination_dv),
                 (3, 'Orbital Period', orbital_period),
                 (4, 'Two-body Calculation', two_body),
                 (5, 'Rocket Efficiency', rocket_equation),
                 (6, 'Orbital Eccentricity', orbital_eccentricity),
                 (7, 'Hyperbolic Excess Velocity', hyperbolic_velocity),
                 (8, 'Create a Custom Object', create_body)]
    print('Select a function:')
    for function in functions:
        print(f'{function[0]}: {function[1]}')
    print("For custom arguments, end input with 'a'. To view documentation, end input with 'd'.")

    prompt = input()
    kwargs = ''
    if prompt[-1] == 'a':
        kwargs = input('Enter argument:\n')
        kwargs = kwargs.replace('true', 'True') # Lazy Caps
        kwargs = kwargs.replace('false', 'False')

    for function in functions:
        if int(prompt[0]) == function[0]:
            print()
            print(f'{prompt}: {function[1]}\n')

            if prompt[-1] == 'd':
                print(function[2].__doc__)
                return True

            if kwargs:
                try:
                    return_value = eval('function[2] (ui=True, ' + kwargs + ')') # dont worry about the eval man
                except SyntaxError:
                    print('SyntaxError: Enter arguments correctly.')
                    print('Tip: Enter quantity with units as string.')
                    return False
            else:
                return_value = function[2]()
    return return_value


if __name__ == '__main__':
    RETURN_VALUE = True
    while RETURN_VALUE:
        RETURN_VALUE = main()
        if RETURN_VALUE:
            input('Press Enter to continue...\n')
