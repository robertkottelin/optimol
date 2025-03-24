import os
from celery import Celery

def make_celery(app=None):
    """Create and configure Celery instance with Flask app context if provided"""
    # Redis connection string with fallback for local development
    redis_url = os.environ.get('REDIS_URL', 'redis://localhost:6379/0')
    
    celery = Celery(
        'optimol_tasks',
        broker=redis_url,
        backend=redis_url,
        include=['tasks']  # Import tasks module where tasks will be defined
    )
    
    # Configure Celery
    celery.conf.update(
        result_expires=3600,  # Results expire after 1 hour
        task_acks_late=True,  # Tasks acknowledged after execution (better error handling)
        task_reject_on_worker_lost=True,  # Allow tasks to be re-queued if worker is lost
        worker_prefetch_multiplier=1,  # Process one task at a time per worker
        task_serializer='json',
        accept_content=['json'],
        result_serializer='json',
        task_routes={
            'tasks.optimize_classical_task': {'queue': 'classical'},
            'tasks.optimize_quantum_task': {'queue': 'quantum'},
            'tasks.optimize_classical_combined_task': {'queue': 'classical'},
            'tasks.optimize_quantum_combined_task': {'queue': 'quantum'},
        }
    )
    
    # Configure Flask integration if app is provided
    if app:
        class FlaskTask(celery.Task):
            def __call__(self, *args, **kwargs):
                with app.app_context():
                    return self.run(*args, **kwargs)
                    
        celery.Task = FlaskTask
    
    return celery