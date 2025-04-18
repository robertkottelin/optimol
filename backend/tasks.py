from celery_config import make_celery
from extensions import db
from sqlalchemy.orm import sessionmaker
import json
import logging
from datetime import datetime
from app import app

import os
import sys
# Add application directory to Python path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from celery_config import make_celery

# Configure logging
logger = logging.getLogger(__name__)

# Initialize Celery
celery = make_celery()

# This will store references to our optimization functions when imported
optimization_funcs = {}

# Task states
TASK_STATUS = {
    'PENDING': 'pending',
    'STARTED': 'running',
    'SUCCESS': 'completed',
    'FAILURE': 'failed',
    'REVOKED': 'cancelled'
}

# Model imports will be done at runtime to avoid circular imports
def import_dependencies():
    """Import dependencies only when needed to avoid circular imports"""
    global optimization_funcs
    
    # Only import if not already imported
    if not optimization_funcs:
        from opti import (
            optimize_classical,
            optimize_quantum,
            optimize_classical_combined,
            optimize_quantum_combined,
            Optimization,
            User,
            ITERATION_LIMITS
        )
        
        # Store function references
        optimization_funcs = {
            'optimize_classical': optimize_classical,
            'optimize_quantum': optimize_quantum,
            'optimize_classical_combined': optimize_classical_combined,
            'optimize_quantum_combined': optimize_quantum_combined,
            'Optimization': Optimization,
            'User': User,
            'ITERATION_LIMITS': ITERATION_LIMITS
        }
        
    return optimization_funcs


@celery.task(bind=True, name='tasks.optimize_classical_task')
def optimize_classical_task(self, atoms, params=None, user_id=None):
    """Celery task for classical optimization"""
    import_dependencies()
    
    logger.info(f"Starting classical optimization task {self.request.id} for user {user_id}")
    
    try:
        with app.app_context():
            # Execute optimization
            result = optimization_funcs['optimize_classical'](atoms, params)
            
            # Store result in database if user_id is provided
            if user_id is not None:
                store_optimization_result(
                    task_id=self.request.id,
                    user_id=user_id, 
                    optimization_type='classical',
                    parameters=json.dumps({
                        'optimization_params': params,
                        'interaction_mode': False,
                        'molecule1_atom_count': len(atoms) if atoms else 0,
                        'molecule2_atom_count': 0
                    }),
                    result=json.dumps(result)
                )
            
            return result
    except Exception as e:
        logger.error(f"Classical optimization task {self.request.id} failed: {str(e)}", exc_info=True)
        raise


@celery.task(bind=True, name='tasks.optimize_quantum_task')
def optimize_quantum_task(self, atoms, params=None, user_id=None):
    """Celery task for quantum optimization"""
    import_dependencies()
    
    logger.info(f"Starting quantum optimization task {self.request.id} for user {user_id}")
    
    try:
        with app.app_context():
            # Execute optimization
            result = optimization_funcs['optimize_quantum'](atoms, params)
            
            # Store result in database if user_id is provided
            if user_id is not None:
                store_optimization_result(
                    task_id=self.request.id,
                    user_id=user_id,
                    optimization_type='quantum',
                    parameters=json.dumps({
                        'optimization_params': params,
                        'interaction_mode': False,
                        'molecule1_atom_count': len(atoms) if atoms else 0,
                        'molecule2_atom_count': 0
                    }),
                    result=json.dumps(result)
                )
            
            return result
    except Exception as e:
        logger.error(f"Quantum optimization task {self.request.id} failed: {str(e)}", exc_info=True)
        raise


@celery.task(bind=True, name='tasks.optimize_classical_combined_task')
def optimize_classical_combined_task(self, molecule1_atoms, molecule2_atoms, params=None, user_id=None, optimize_molecule1=True, optimize_molecule2=True):
    """Celery task for classical combined optimization"""
    import_dependencies()
    
    logger.info(f"Starting classical combined optimization task {self.request.id} for user {user_id}, optimize_molecule1={optimize_molecule1}, optimize_molecule2={optimize_molecule2}")
    
    try:
        with app.app_context():
            # Execute optimization with the additional parameters
            result = optimization_funcs['optimize_classical_combined'](
                molecule1_atoms, 
                molecule2_atoms, 
                params, 
                optimize_molecule1, 
                optimize_molecule2
            )
            
            # Store result in database if user_id is provided
            if user_id is not None:
                store_optimization_result(
                    task_id=self.request.id,
                    user_id=user_id,
                    optimization_type='classical',
                    parameters=json.dumps({
                        'optimization_params': params,
                        'interaction_mode': True,
                        'molecule1_atom_count': len(molecule1_atoms) if molecule1_atoms else 0,
                        'molecule2_atom_count': len(molecule2_atoms) if molecule2_atoms else 0,
                        'optimize_molecule1': optimize_molecule1,
                        'optimize_molecule2': optimize_molecule2
                    }),
                    result=json.dumps(result)
                )
            
            return result
    except Exception as e:
        logger.error(f"Classical combined optimization task {self.request.id} failed: {str(e)}", exc_info=True)
        raise


@celery.task(bind=True, name='tasks.optimize_quantum_combined_task')
def optimize_quantum_combined_task(self, molecule1_atoms, molecule2_atoms, params=None, user_id=None, optimize_molecule1=True, optimize_molecule2=True):
    """Celery task for quantum combined optimization"""
    import_dependencies()
    
    logger.info(f"Starting quantum combined optimization task {self.request.id} for user {user_id}, optimize_molecule1={optimize_molecule1}, optimize_molecule2={optimize_molecule2}")
    
    try:
        with app.app_context():
            # Execute optimization with selective molecule optimization parameters
            result = optimization_funcs['optimize_quantum_combined'](
                molecule1_atoms, 
                molecule2_atoms, 
                params, 
                optimize_molecule1, 
                optimize_molecule2
            )
            
            # Store result in database if user_id is provided
            if user_id is not None:
                store_optimization_result(
                    task_id=self.request.id,
                    user_id=user_id,
                    optimization_type='quantum',
                    parameters=json.dumps({
                        'optimization_params': params,
                        'interaction_mode': True,
                        'molecule1_atom_count': len(molecule1_atoms) if molecule1_atoms else 0,
                        'molecule2_atom_count': len(molecule2_atoms) if molecule2_atoms else 0,
                        'optimize_molecule1': optimize_molecule1,
                        'optimize_molecule2': optimize_molecule2
                    }),
                    result=json.dumps(result)
                )
            
            return result
    except Exception as e:
        logger.error(f"Quantum combined optimization task {self.request.id} failed: {str(e)}", exc_info=True)
        raise

def store_optimization_result(task_id, user_id, optimization_type, parameters, result):
    """Store optimization result in database"""
    import_dependencies()
    
    # Create new SQLAlchemy session
    Session = sessionmaker(bind=db.engine)
    session = Session()
    
    try:
        with app.app_context():
            # Find existing record by task_id
            optimization = session.query(optimization_funcs['Optimization']).filter_by(id=task_id).first()
            
            if optimization:
                # Update existing record with results
                optimization.result = result
                logger.info(f"Updated existing optimization result for task {task_id}")
            else:
                # Fallback: create new record if not found
                optimization = optimization_funcs['Optimization'](
                    id=task_id,
                    user_id=user_id,
                    optimization_type=optimization_type,
                    parameters=parameters,
                    result=result,
                    created_at=datetime.now()
                )
                session.add(optimization)
                logger.info(f"Created new optimization result for task {task_id}")
                
            session.commit()
    except Exception as e:
        logger.error(f"Failed to store optimization result for task {task_id}: {str(e)}", exc_info=True)
        session.rollback()
        raise
    finally:
        session.close()