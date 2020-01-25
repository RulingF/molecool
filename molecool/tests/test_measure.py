"""
Unit tests for measure module.
"""

import numpy as np
import pytest
import molecool

def test_calculate_distance():
    """Test that calculate_distance function calculates what we expect."""

    r1 = np.array([0,0,0])
    r2 = np.array([0,0.1,0])

    expected_distance = 0.1

    calculated_distance = molecool.calculate_distance(r1,r2)

    assert calculated_distance == expected_distance

def test_calculate_distance_typeerror():
    """Test type error"""

    r1 = [0,0,0]
    r2 = [1,0,0]
#    with pytest.raises(TypeError,'np.ndarray'):
    with pytest.raises(TypeError):
        calculate_distance = molecool.calculate_distance(r1,r2)

def test_calculate_angle():
    """Test that calculate_angle function calcualtes what we expect."""

    r1 = np.array([0,0,-1])
    r2 = np.array([0,0,0])
    r3 = np.array([1,0,0])

    expected_angle = 90.

    calculated_angle = molecool.calculate_angle(r1,r2,r3,degrees = True)

    assert calculated_angle == expected_angle

#@pytest.mark.parametrize("variable_name1, variable_name2, variable_name3, ..., variable_nameN, expected_answer", [
#    (variable_value1, variable_value2, variable_value3, ..., variable_valueN, expected_answer),
#    (variable_value1, variable_value2, variable_value3, ..., variable_valueN, expected_answer),
#]) 
#def test_name(variable_value1, variable_value2, variable_value3, ..., variable_valueN, expected_answer):

@pytest.mark.parametrize("p1,p2,p3, expected_angle", [
    (np.array([np.sqrt(2)/2,np.sqrt(2)/2, 0]), np.array([0,0,0]), np.array([1,0,0]),45),
    (np.array([0,0,-1]),np.array([0,1,0]),np.array([1,0,0]),60)
])
def test_calculate_angle_many(p1,p2,p3,expected_angle):
    """Test many angle values"""

    calculated_angle = molecool.calculate_angle(p1,p2,p3,degrees = True)

    assert expected_angle == pytest.approx(calculated_angle)

     
